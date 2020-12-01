/*
* Copyright 2018-2020 Redis Labs Ltd. and Contributors
*
* This file is available under the Redis Labs Source Available License Agreement
*/

#include "execution_ctx.h"
#include "../query_ctx.h"
#include "../execution_plan/execution_plan_clone.h"

#include "lucata_state.h"

#include <pcre.h>

LucataState *g_lucataState = NULL;

static ExecutionType _GetExecutionTypeFromAST(AST *ast) {
	const cypher_astnode_type_t root_type = cypher_astnode_type(ast->root);
	if(root_type == CYPHER_AST_QUERY) return EXECUTION_TYPE_QUERY;
	if(root_type == CYPHER_AST_CREATE_NODE_PROPS_INDEX) return EXECUTION_TYPE_INDEX_CREATE;
	if(root_type == CYPHER_AST_DROP_NODE_PROPS_INDEX) return EXECUTION_TYPE_INDEX_DROP;
	assert(false && "Uknown execution type");
}

static ExecutionCtx *_ExecutionCtx_New(AST *ast, ExecutionPlan *plan, ExecutionType exec_type) {
	ExecutionCtx *exec_ctx = rm_calloc(1, sizeof(ExecutionCtx));
	exec_ctx->ast = ast;
	exec_ctx->plan = plan;
	exec_ctx->exec_type = exec_type;
	return exec_ctx;
}

static ExecutionCtx _ExecutionCtx_Clone(const ExecutionCtx orig) {
	ExecutionCtx ctx = {0};
	ctx.ast = AST_ShallowCopy(orig.ast);
	ctx.plan = ExecutionPlan_Clone(orig.plan);
	ctx.exec_type = orig.exec_type;
	return ctx;
}

static ExecutionCtx *_ExecutionCtx_FromCache(Cache *cache, const char *query_string,
											 cypher_parse_result_t *params_parse_result) {
	// Check the cache to see if we already have a cached context for this query.
	ExecutionCtx *cached_exec_ctx = Cache_GetValue(cache, query_string);
	if(cached_exec_ctx) {
		// Cache hit - Clone the execution context. Set the execution type for query execution and indicate a cache hit.
		// Set AST as it is retrived from cache and it is required for execution plan clone.
		QueryCtx_SetAST(cached_exec_ctx->ast);
	}
	return cached_exec_ctx;
}

static AST *_ExecutionCtx_ParseAST(const char *query_string,
								   cypher_parse_result_t *params_parse_result) {
	cypher_parse_result_t *query_parse_result = parse_query(query_string);
	// If no output from the parser, the query is not valid.
	if(!query_parse_result) {
		parse_result_free(params_parse_result);
		return NULL;
	}

	// Prepare the constructed AST.
	AST *ast = AST_Build(query_parse_result);
	// Set parameters parse result in the execution ast.
	AST_SetParamsParseResult(ast, params_parse_result);
	return ast;
}

ExecutionCtx ExecutionCtx_FromQuery(const char *query) {
    const char *regex = "MATCH \\(s:([^)]+)\\)-\\[\\*([^)]+)\\]->\\(t\\) WHERE s.id=([^)]+) RETURN count\\(t\\)";
    //printf("%s query: %s[%lu] regex: %s [%lu]\n", __func__, query, strlen(query), regex, strlen(regex));
    pcre *reCompiled = NULL;
    const char *pcreErrorStr;
    int pcreErrorOffset;
    reCompiled = pcre_compile(regex, 0, &pcreErrorStr, &pcreErrorOffset, NULL);
    if (reCompiled) {
//        printf("Succesfully compiled regex into pattern buffer\n");
    } else {
        printf("Failed to compile regex into pattern buffer\n");
    }

    pcre_extra *pcreExtra = pcre_study(reCompiled, 0, &pcreErrorStr);
    if(pcreErrorStr != NULL) {
        printf("ERROR: Could not study '%s': %s\n", regex, pcreErrorStr);
    }

	int subStrVec[30];
	int pcreExecRet = pcre_exec(reCompiled,
                            pcreExtra,
                            query, 
                            strlen(query),  // length of string
                            0,                      // Start looking at this point
                            0,                      // OPTIONS
                            subStrVec,
                            30);
    if (pcreExecRet < 0) {
  //      printf("NO MATCH FOUND %d\n", pcreExecRet);

    } else {
        g_lucataState = (LucataState *)malloc(sizeof(LucataState));
        g_lucataState->m_runningKhop = true;
        g_lucataState->m_khopResultsAvailable = false;
        g_lucataState->m_khopResultsObtained = false;
  //      printf("MATCH FOUND! %d\n", pcreExecRet);
        // 0 is whole, 1 is graph id, 2 is number of hops, 3 is seed id
        for(int j=0; j < pcreExecRet; j++) {
			const char *psubStrMatchStr;

        	pcre_get_substring(query, subStrVec, pcreExecRet, j, &(psubStrMatchStr));
   //     	printf("Match(%2d/%2d): (%2d,%2d): '%s'\n", j, pcreExecRet-1, subStrVec[j*2], subStrVec[j*2+1], psubStrMatchStr);
      	}
    }

	// Have an invalid ctx for errors.
	ExecutionCtx invalid_ctx = {.ast = NULL, .plan = NULL, .cached = false, .exec_type = EXECUTION_TYPE_INVALID};
	const char *query_string;
	// Parse and validate parameters only. Extract query string.
	cypher_parse_result_t *params_parse_result = parse_params(query, &query_string);
	// Return invalid execution context if there isn't a parser result.
	if(params_parse_result == NULL) return invalid_ctx;

	GraphContext *gc = QueryCtx_GetGraphCtx();
	Cache *cache = GraphContext_GetCache(gc);
	// Check the cache to see if we already have a cached context for this query.
	ExecutionCtx *cached_exec_ctx = _ExecutionCtx_FromCache(cache, query_string, params_parse_result);
	if(cached_exec_ctx) {
		ExecutionCtx ctx = _ExecutionCtx_Clone(*cached_exec_ctx);
		// Set parameters parse result in the execution ast.
		AST_SetParamsParseResult(ctx.ast, params_parse_result);
		ctx.cached = true;
		return ctx;
	}

	// No cached execution plan, try to parse the query.
	AST *ast = _ExecutionCtx_ParseAST(query_string, params_parse_result);
	// Invalid query, return invalid execution context.
	if(!ast) return invalid_ctx;

	ExecutionPlan *plan = NULL;
	ExecutionType exec_type = _GetExecutionTypeFromAST(ast);
	// In case of valid query, create execution plan, and cache it and the AST.
	if(exec_type == EXECUTION_TYPE_QUERY) {
		plan = NewExecutionPlan();
		// Created new valid execution context.
		ExecutionCtx *exec_ctx_to_cache = _ExecutionCtx_New(ast, plan, exec_type);
		// Cache execution context.
		Cache_SetValue(cache, query_string, exec_ctx_to_cache);
		// Clone execution plan and ast that will be used in the current execution.
		plan = ExecutionPlan_Clone(plan);
		ast = AST_ShallowCopy(ast);
	}
	ExecutionCtx ctx = {.ast = ast, .plan = plan, .exec_type = exec_type, .cached = false};
	return ctx;
}

void ExecutionCtx_Free(ExecutionCtx *ctx) {
	if(!ctx) return;
	if(ctx->ast) AST_Free(ctx->ast);
	if(ctx->plan) ExecutionPlan_Free(ctx->plan);
	rm_free(ctx);
}

