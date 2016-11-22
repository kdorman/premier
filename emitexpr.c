#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "emitexpr.h"
#include "message.h"
#include "kmer.h"

static const int emitexpr_token2op_cmp[] = {
	OP_CMP_LT, OP_CMP_LE, OP_CMP_GT, OP_CMP_GE
};

const function_t emitexpr_func_table[4] = {
	{"all", 0, T_BOOLS, T_INDICATOR},
	{"any", 1, T_BOOLS, T_INDICATOR},
	{"count", 2, T_BOOLS, T_COUNTER},
	{0, 0, 0, 0}
};

const char *emitexpr_token_names[] = {
	"T_OP_L_SUBSTR", "T_OP_R_SUBSTR", "T_OP_FUNC_NAME", 
	"T_OP_L_PAREN", "T_OP_R_PAREN", "T_OP_EQ", "T_OP_LT", 
	"T_OP_LE", "T_OP_GT", "T_OP_GE", "T_OP_NE",
	"T_OP_PLUS", "T_OP_MINUS",
	"T_OP_SLICE", "T_OP_COMMA", "T_OP_SUBSTR_OF", 
	"T_OP_FUNCTION", 
	"T_STR_LITERAL", "T_NUM_LITERAL", 
	"T_SLICE", "T_XT", "T_ST", "T_YT", "T_BASES", 
	"T_KMER_LEN","T_BOOLS", "T_INDICATOR", "T_COUNTER",  
};

const int emitexpr_precedence_table[] = {
	/* T_OP_L_SUBSTR, T_OP_R_SUBSTR */
	0, 1,
	/* T_OP_L_FUNC_PAREN, T_OP_L_PAREN, T_OP_R_PAREN */
	0, 0, 1,
	/* T_OP_EQ, T_OP_LT, T_OP_LE, T_OP_GT, T_OP_GE, T_OP_NE */
	4, 4, 4, 4, 4, 5,
	/* T_OP_PLUS, T_OP_MINUS */
	5, 5,
	/* T_OP_SLICE */
	3,
	/* T_OP_COMMA */
	2,
	/* T_OP_SUBSTR_OF */ 
	6,
	/* T_OP_FUNCTION */
	6
};

const int emitexpr_assoc_table[] = {
	ASSOC_NONE,	ASSOC_NONE, 	
	ASSOC_NONE,	ASSOC_NONE, ASSOC_NONE,
	/* T_OP_EQ, T_OP_LT, T_OP_LE, T_OP_GT, T_OP_GE, T_OP_NE */
	ASSOC_BINARY, ASSOC_BINARY, ASSOC_BINARY, ASSOC_BINARY, ASSOC_BINARY, ASSOC_BINARY,
	/* T_OP_PLUS, T_OP_MINUS */
	ASSOC_BINARY, ASSOC_BINARY,
	/* T_OP_SLICE */
	ASSOC_BINARY,
	/* T_OP_COMMA */
	ASSOC_BINARY,
	/* T_OP_SUBSTR_OF */ 
	ASSOC_UNITARY_RIGHT,
	/* T_OP_FUNCTION */
	ASSOC_UNITARY_RIGHT
};

static int emitexpr_reduce(stack_token_t *st, stack_ast_t *sa, 
		int max_prec, int token_prec, size_t *pradix, size_t *pinit_radix, 
		int kmer_len, int qmax);
static int emitexpr_lexer(char *pc, char **pe);
static void emitexpr_build_ast(token_t *token, stack_ast_t *sa, int discard);
static void emitexpr_ast_traverse(ast_node_t *node, int idx,
		void (*callback) (ast_node_t *node, int idx, void *data), 
		void *data, int to_free);
static void emitexpr_cb_compile(ast_node_t *node, int idx, void *data);
static void emitexpr_push_back_op(oplist_t *oplst, uint8_t op_code, 
		uint8_t src1, uint8_t src2, uint8_t dst, void *param);

/** Traverse a built AST using user designated callback function */
static void emitexpr_ast_traverse(ast_node_t *node, int idx,
		void (*callback) (ast_node_t *node, int idx, void *data), 
		void *data, int to_free)
{
	if (node->lhs) {
		emitexpr_ast_traverse(node->lhs, idx*2+1, callback, data, to_free);
	}

	if (node->rhs) {
		emitexpr_ast_traverse(node->rhs, idx*2+2, callback, data, to_free);
	}

	callback(node, idx, data);

	if (to_free) free(node);
}

static inline void emitexpr_push_back_op(oplist_t *oplst, uint8_t op_code, 
		uint8_t src1, uint8_t src2, uint8_t dst, void *param)
{
	op_t new_op = (op_t) {op_code, src1, src2, dst, param};
	oplst->ops[oplst->n_ops++] = new_op;
}

static inline int _find_avail_reg(uint64_t *reg_usage) {
	register uint64_t _reg = *reg_usage;
	for (int i = 0; i < N_REGISTERS; i++) {
		if (((_reg >> i) & 1) == 0) {
			*reg_usage = _reg | (1 << i);
			return i;
		}
	}
	return 0;
}

static inline uint64_t _eval_yt(char *yt, int yt_len, int qmax)
{
	uint64_t ynum = 0;
	uint64_t radix = 1;
	for (int i = 0; i < yt_len; i++) {
		ynum += radix * yt[i];
		radix *= qmax;
	}

	return ynum;
}

#define _reset_reg(r, i) ((r) &= ~(1UL << i))

static void emitexpr_cb_compile(ast_node_t *node, int idx, void *data)
{
	oplist_t *oplst = data;
	stack_register_t *sr = oplst->sr;

	static uint64_t reg_usage = 0;
	int src1, src2, dst;

	token_t *token = &node->token;

	switch (token->token_type) {
		case T_NUM_LITERAL:
			dst = _find_avail_reg(&reg_usage);
			stack_push(*sr, dst);

			INFO("MOVE %lu TO $%d\n", token->immediate, dst);
			emitexpr_push_back_op(oplst, OP_MOVE,
						0xFF, 0xFF, dst, (void *) token->immediate);

			break;
		case T_ST:
		case T_XT:
		case T_YT:
			dst = _find_avail_reg(&reg_usage);
			stack_push(*sr, dst);

			INFO("LOAD %s TO $%d\n", 
					token->token_type == T_ST ? "s_t" : (
						token->token_type == T_XT ? "x_t" : "y_t"), dst);

			do {
				emitexpr_push_back_op(oplst, 
						token->token_type == T_ST ? OP_LOAD_ST : (
							token->token_type == T_XT ? OP_LOAD_XT : 
							OP_LOAD_YT),
						0xFF, 0xFF, dst, NULL);

				int ub = (int) (token->immediate >> 32);
				int lb = (int) (token->immediate & UINT32_MAX);
				uint64_t ctx_len =  ub - lb + 1;

				if (token->token_as == T_BASES) {
					emitexpr_push_back_op(oplst, 
							OP_SUBSTRQ,
							dst, 0xFF, dst, (void *) token->immediate);

					INFO("SUBSTRQ %d %d ON $%d TO $%d\n", 
							lb, ub, dst, dst);
				}
				else {
					/*
					emitexpr_push_back_op(oplst,
							OP_TYPE, 
							0xFF, 0xFF, dst, 
							(void *) T_YT);
					*/
					emitexpr_push_back_op(oplst,
							OP_ADD,
							dst, 0xFF, dst, 
							(void *) lb); 
					emitexpr_push_back_op(oplst,
							OP_SET_YT_LEN,
							0xFF, 0xFF, 0xFF, (void *) ctx_len);

					/*
					INFO("TYPE %d ON $%d\n",
							T_YT, dst);
					*/
					INFO("ADD %d ON $%d TO $%d\n", 
							lb, dst, dst);
					INFO("SET_YT_LEN %d\n", ctx_len);
				}

				emitexpr_push_back_op(oplst, 
						OP_PUSH,
						0xFF, 0xFF, 0xFF, (void *) ctx_len);

				INFO("PUSH %lu\n", ctx_len);

			} while (0);

			break;

		case T_OP_FUNCTION:
			src1 = stack_pop(*sr);
			dst = src1;
			stack_push(*sr, dst);

			emitexpr_push_back_op(oplst,
					OP_CALL, 
					src1, 0xFF, dst, (void *) token->immediate);

			emitexpr_push_back_op(oplst,
					OP_TYPE, 
					0xFF, 0xFF, dst, 
					(void *) emitexpr_func_table[
						node->token.immediate].return_type);

			INFO("TYPE %d ON $%d\n",
					emitexpr_func_table[
						node->token.immediate].return_type,
						dst);

			INFO("CALL %s ON $%d TO $%d\n", 
					emitexpr_func_table[node->token.immediate].name,
					src1, dst);
			break;

		case T_OP_EQ:
		case T_OP_NE:
			src2 = stack_pop(*sr);
			src1 = stack_pop(*sr);

			emitexpr_push_back_op(oplst,
					token->token_type == T_OP_EQ ?  OP_CMP_EQ : OP_CMP_NE,
					src1, src2, src1, NULL);

			INFO("%s $%d WITH $%d TO $%d\n",
					token->token_type == T_OP_EQ ? "CMPEQ" : "CMPNE",
					src1, src2, src1);

			stack_push(*sr, src1);

			emitexpr_push_back_op(oplst,
					OP_POP,
					0xFF, 0xFF, src2, NULL);

			INFO("POP TO $%d\n", src2);

			_reset_reg(reg_usage, src2);

			break;

		case T_OP_LT:
		case T_OP_LE:
		case T_OP_GT:
		case T_OP_GE:
			src2 = stack_pop(*sr);
			src1 = stack_pop(*sr);

			do {
				int op_code = emitexpr_token2op_cmp[
					token->token_type - T_OP_LT];
				/* the left operand HAS TO BE T_YT, because if both operands
				 * are immediates, they should have been reduced already. */	
				emitexpr_push_back_op(oplst,
						op_code,	
						src1, src2, src1, NULL);
			} while (0);

			INFO("CMP%s $%d WITH $%d TO $%d\n",
					emitexpr_token_names[token->token_type] + 5,
					src1, src2, src1);

			stack_push(*sr, src1);
			_reset_reg(reg_usage, src2);

			/* for correctly calling function */
			INFO("PEEK STACK\n");
			emitexpr_push_back_op(oplst,
					OP_PEEK,	
					0xFF, 0xFF, 0xFF, NULL);

			break;

		case T_OP_COMMA:
			dst = _find_avail_reg(&reg_usage);
			INFO("POP TO $%d\n", dst);
			INFO("RADIX AT $%d\n", dst);

			emitexpr_push_back_op(oplst,
					OP_POP,
					0xFF, 0xFF, dst, NULL);
			emitexpr_push_back_op(oplst,
					OP_RADIX_SET, 
					dst, 0xFF, 0xFF, NULL);

			src2 = stack_pop(*sr);
			src1 = stack_pop(*sr);

			emitexpr_push_back_op(oplst,
					OP_REDUCE,
					src1, src2, src1, NULL);

			INFO("REDUCE $%d WITH $%d TO $%d\n",
					src1, src2, src1);

			stack_push(*sr, src1);

			_reset_reg(reg_usage, src2);
			_reset_reg(reg_usage, dst);

			break;
	}
}

static void emitexpr_build_ast(token_t *token, stack_ast_t *sa, int discard)
{
	int assoc = ASSOC_OF(token->token_type);
	ast_node_t *new_ast_node = malloc(sizeof(*new_ast_node));
	new_ast_node->token = *token;

	if (!discard) {
		if (assoc == ASSOC_BINARY) {
			new_ast_node->rhs = stack_pop(*sa);
			new_ast_node->lhs = stack_pop(*sa);
		}
		else if (assoc == ASSOC_UNITARY_RIGHT) {
			new_ast_node->lhs = stack_pop(*sa);
			new_ast_node->rhs = NULL;
		}
	}
	else {
		new_ast_node->lhs = NULL;
		new_ast_node->rhs = NULL;

		if (assoc == ASSOC_BINARY) {
			free(stack_pop(*sa));
			free(stack_pop(*sa));
		}
		else if (assoc == ASSOC_UNITARY_RIGHT) {
			free(stack_pop(*sa));
		}
	}

	stack_push(*sa, new_ast_node);
}

static int emitexpr_reduce(stack_token_t *st, stack_ast_t *sa, 
		int max_prec, int token_prec, size_t *pradix, size_t *pinit_radix, 
		int kmer_len, int qmax)
{
	int new_max_prec = max_prec;
	function_t const *f;

	if (stack_empty(*st)) {
		return -1;
	}

	while (new_max_prec > token_prec) {
		token_t *t_op = stack_top(*st);
		for (; !IS_OP(t_op->token_type) || t_op->precedence <= token_prec;
				t_op--) {
			if (t_op == stack_bottom(*st)) break; 
		}

		switch (t_op->token_type) {
			/* since we don't have any variables, all arithmatic operations
			 * involve with numeric constants, so simply reduce the expression
			 * to another constant */
			case T_OP_LT:
			case T_OP_LE:
			case T_OP_GT:
			case T_OP_GE:
				/* return T_BOOLS if any of the operand is T_YT, 
				 * or T_NUM_LITERAL if both operands are T_NUM_LITERAL */
				if (t_op == stack_bottom(*st) || t_op == stack_top(*st)) {
					INFO("Insufficient operands!\n");
					exit(1);
				}

				if ((t_op-1)->token_type == T_NUM_LITERAL && 
						(t_op+1)->token_type == T_NUM_LITERAL) {

					uint64_t imm_a = (t_op-1)->immediate;
					uint64_t imm_b = (t_op+1)->immediate;

					(t_op-1)->immediate = t_op->token_type == T_OP_LT ? 
						imm_a < imm_b : (t_op->token_type == T_OP_LE ? 
								imm_a <= imm_b : (t_op->token_type == T_OP_GT ?
									imm_a > imm_b : imm_a >= imm_b));

					emitexpr_build_ast(t_op, sa, 1);

					/* pop out rhs and op */
					stack_pop(*st);
					stack_pop(*st);
				}
				else if ((t_op-1)->token_type == T_YT && 
						(t_op+1)->token_type == T_NUM_LITERAL) {
					
					(t_op-1)->token_type = T_BOOLS;

					INFO("[CMP] Compare y_t and %lu\n", (t_op+1)->immediate);
					emitexpr_build_ast(t_op, sa, 0);

					stack_pop(*st);
					stack_pop(*st);
				}

				break;

			case T_OP_PLUS:
			case T_OP_MINUS:
				if (t_op == stack_bottom(*st) ||
						t_op == stack_top(*st) ||
						(t_op-1)->token_type != T_NUM_LITERAL ||
						(t_op+1)->token_type != T_NUM_LITERAL) {
					INFO("Invalid numerical values\n");
					exit(1);
				}

				if (t_op->token_type == T_OP_MINUS) {
					int c = (int) (t_op-1)->immediate - (int) (t_op+1)->immediate;

					INFO("[MINUS] %d - %d = %d\n",
							(int) (t_op-1)->immediate, (int) (t_op+1)->immediate, c);

					(t_op-1)->immediate = c;

					stack_pop(*sa);
					stack_pop(*sa);
					emitexpr_build_ast(t_op-1, sa, 1);

					/* pop out rhs and op */
					stack_pop(*st);
					stack_pop(*st);
				}

				break;
			case T_OP_SLICE:
				if (t_op == stack_bottom(*st) ||
						t_op == stack_top(*st) ||
						(t_op-1)->token_type != T_NUM_LITERAL ||
						(t_op+1)->token_type != T_NUM_LITERAL) {
					INFO("Invalid slicing lower and upper bound\n");
					exit(1);
				}

				do {
					INFO("[SLICE] %lu..%lu\n",
							(t_op-1)->immediate, (t_op+1)->immediate);
					uint64_t slice_bound = ((t_op+1)->immediate << 32) | 
						(t_op-1)->immediate;

					emitexpr_build_ast(t_op, sa, 1);

					stack_pop(*st);
					stack_pop(*st);

					token_t *tk_slice = stack_top(*st);
					tk_slice->token_type = T_SLICE;
					tk_slice->immediate = slice_bound;

					(*stack_top(*sa))->token = *(tk_slice);

				} while (0);

				break;

			case T_OP_EQ:
			case T_OP_NE:
				if (t_op == stack_bottom(*st) ||
						t_op == stack_top(*st) ||
						(t_op-1)->token_as != T_BASES ||
						(t_op+1)->token_as != T_BASES) {
					INFO("Can only compare contexts\n");
					exit(1);
				}

				if ((t_op-1)->immediate != (t_op+1)->immediate) {
					INFO("Context ranges don't match!\n");
					exit(1);
				}
				
				INFO("[CMP] Compare %s and %s\n",
						(t_op-1)->token_type == T_ST ? "s_t" : "x_t",
						(t_op+1)->token_type == T_ST ? "s_t" : "x_t");

				(t_op-1)->token_type = T_BOOLS;

				emitexpr_build_ast(t_op, sa, 0);

				stack_pop(*st);
				stack_pop(*st);

				break;

			case T_OP_COMMA:
				if (t_op == stack_bottom(*st) || 
						t_op == stack_top(*st)) {
					INFO("Need both lhs and rhs for comma.\n");
					exit(1);
				}

				if ((t_op-1)->token_as != T_BASES && 
						(t_op-1)->token_type != T_COUNTER &&
						(t_op-1)->token_type != T_INDICATOR &&
						(t_op-1)->token_type != T_YT) {
					INFO("The left operand should be typed as bases, counter,"
							" indicator or quality scores [observed: %s]\n",
							emitexpr_token_names[(t_op-1)->token_type]);
					exit(1);
				}

				if ((t_op+1)->token_as != T_BASES && 
						(t_op+1)->token_type != T_COUNTER &&
						(t_op+1)->token_type != T_INDICATOR &&
						(t_op+1)->token_type != T_YT) {
					INFO("The left operand should be typed as bases, counter,"
							" indicator or quality scores [observed: %s]\n",
							emitexpr_token_names[(t_op+1)->token_type]);
					exit(1);
				}

				INFO("[COMMA] Computing context index for %s and %s\n",
						(t_op-1)->token_as == T_BASES ? "BASES" : 
						((t_op-1)->token_type == T_COUNTER ? "COUNTER" : 
						 ((t_op-1)->token_type == T_INDICATOR ? "INDICATOR" :
						  "QSCORES")),
						(t_op+1)->token_as == T_BASES ? "BASES" : 
						((t_op+1)->token_type == T_COUNTER ? "COUNTER" : 
						 ((t_op+1)->token_type == T_INDICATOR ? "INDICATOR" :
						  "QSCORES")));
				do {
					token_t *tk_rhs = t_op + 1;
					int ctx_len = (tk_rhs->immediate >> 32) - 
								(tk_rhs->immediate & UINT32_MAX) + 1;

					if (tk_rhs->token_type == T_YT && ctx_len > 1) {
						INFO("WARNING: Using y_t in the conditioning can "
								"quickly explode the number of PMFs! "
								"Consider using indicators or counters!\n");
					}

					if (tk_rhs->token_type == T_ST && 
							(tk_rhs->immediate >> 32) == kmer_len) {
						*pinit_radix = *pradix * (int) pow(4, ctx_len-1);
					}

					/* if an indicator, can have up to two values: 0, 1,
					 * for a counter, can have up to ctx_len+1 values: 
					 *		0, 1, * ..., ctx_len
					 * and finally, for a substring of length ctx_len,
					 * can have 4^{ctx_len} number of possibilities.
					 */
					int rhs_radix = tk_rhs->token_as == T_BASES ? 
						(int) pow(4, ctx_len) : 
						(tk_rhs->token_type == T_COUNTER ? (ctx_len + 1) : 
						 (tk_rhs->token_type == T_INDICATOR ? 2 :
							(int) pow(qmax, ctx_len) 
						 ));

					*pradix *= rhs_radix;
				} while(0);

				emitexpr_build_ast(t_op, sa, 0);

				/* pop comma and the rhs */
				stack_pop(*st);
				stack_pop(*st);

				break;

			case T_OP_R_SUBSTR:
				if (t_op == stack_bottom(*st)) {
					INFO("Unmatched ']'.");
					exit(1);
				} if ((t_op-1) == stack_bottom(*st) ||
						((t_op-1)->token_type != T_SLICE && 
						 (t_op-1)->token_type != T_NUM_LITERAL)) {
					INFO("Substring index has to be a slice or a number.\n");
					exit(1);
				}

				if ((t_op-2)->token_type != T_OP_L_SUBSTR) {
					INFO("Unmatched ']'.\n");
					exit(1);
				}

				do {
					token_t *tk_idx = t_op - 1;
					if ((tk_idx->immediate >> 32) > kmer_len) {
						INFO("Invalid upper bound for substring index.\n");
						exit(1);
					}
				} while (0);

				if ((t_op-1)->token_type == T_NUM_LITERAL) {
					uint64_t idx = (t_op-1)->immediate;
					(t_op-1)->immediate = (idx << 32) | idx;
				}

				(t_op-2)->precedence = PRECEDENCE_OF(T_OP_SUBSTR_OF);
				(t_op-2)->token_type = T_OP_SUBSTR_OF;
				(t_op-2)->immediate = (t_op-1)->immediate;

				INFO("[SUBSTR] Substring from %d to %d\n",
						(int) ((t_op-2)->immediate & UINT32_MAX), 
						(int) ((t_op-2)->immediate >> 32));

				/* immediate reduce within AST : replace T_NUM_LITERAL
				 * or T_SLICE with T_OP_SUBSTR_OF */
				emitexpr_build_ast((t_op-2), sa, 1);

				stack_pop(*st);
				stack_pop(*st);

				return PRECEDENCE_OF(T_OP_SUBSTR_OF);

				break;

			case T_OP_FUNCTION:
				t_op->precedence = 0;
				t_op->token_type = t_op->token_as;

				INFO("[FUNC] Invoking function...\n"); 

				break;

			case T_OP_SUBSTR_OF:
				if (t_op == stack_bottom(*st)) {
					INFO("Dangling substring operator.\n");
					exit(1);
				}

				if ((t_op-1)->token_type != T_ST && 
						(t_op-1)->token_type != T_XT &&
						(t_op-1)->token_type != T_YT) {
					INFO("Invalid subject to get substring of.\n");
					exit(1);
				}

				INFO("[SUBSTR_OF] Substring of %s from %d to %d\n",
						(t_op-1)->token_type == T_ST ? "s_t" : 
						((t_op-1)->token_type == T_XT ? "x_t": "y_t"),
						(int) (t_op->immediate & UINT32_MAX),
						(int) (t_op->immediate >> 32));

				(t_op-1)->immediate = t_op->immediate;
				
				/* immediate reduce within AST */
				stack_pop(*sa);
				(*stack_top(*sa))->token = *(t_op-1);

				stack_pop(*st);

				break;

			case T_OP_R_PAREN:
				if (t_op == stack_bottom(*st)) {
					INFO("Unmatched ')'.\n");
					exit(1);
				}

				if ((t_op-1) == stack_bottom(*st)) {
					INFO("No functional argument.\n");
					exit(1);
				}

				if (((t_op-2)->token_type != T_OP_L_PAREN)) {
					INFO("Unmatched ')'.\n");
					exit(1);
				}

				if ((t_op-2) != stack_bottom(*st) &&
						(t_op-3)->token_type == T_OP_FUNC_NAME) {
					int func_idx = (t_op-3)->immediate;

					f = emitexpr_func_table;
					for (; f->func_index != func_idx; f++);
					if (f->arg_type != (t_op-1)->token_type) {
						INFO("Incompatible functional argument.\n");
						exit(1);
					}

					INFO("[FUNC_PAREN] Checking function invocation `%s`\n",
							f->name);

					token_t *t_func = t_op - 3;
					t_func->precedence = PRECEDENCE_OF(T_OP_FUNCTION);
					t_func->token_type = T_OP_FUNCTION;

					emitexpr_build_ast(t_func, sa, 0);

					t_func->immediate = (t_op-1)->immediate;

					stack_pop(*st);
					stack_pop(*st);
					stack_pop(*st);

					return PRECEDENCE_OF(T_OP_SUBSTR_OF);
				}
				else {
					*(t_op-2) = *(t_op-1);
					stack_pop(*st);
					stack_pop(*st);
				}

				break;
			default:
				return token_prec;
		}

		int updated_flag = 0;
		token_t *tk = stack_top(*st);
		for (int i = st->used; i > 0; i--, tk--) {
			if (IS_OP(tk->token_type)){
				if (tk->precedence <= new_max_prec &&
						tk->precedence >= token_prec) {
					updated_flag = 1;
					new_max_prec = tk->precedence;
					break;
				}
			}
		}

		if (!updated_flag)
			return token_prec;
	}

	return new_max_prec;
}

ast_node_t *emitexpr_parse(char *expr, const char *end_chars, 
		char **endexpr, token_t *reduced_token, size_t *pradix, 
		size_t *pinit_radix, int kmer_len, int qmax)
{
	char *pc = (char *) expr;
	char *pe = NULL;

	int max_precedence = 0;

	token_t redu_tk;

	token_t invalid_tk = {0, T_INVALID, T_INVALID, 0, 0, NULL};
	ast_node_t invalid_ast = {invalid_tk, NULL, NULL};

	stack_token_t expr_stack;
	stack_ast_t ast_stack;

	stack_init(token_t, expr_stack, 64, invalid_tk);
	stack_init(ast_node_t *, ast_stack, 64, &invalid_ast);

	while (*pc != '\0' && strchr(end_chars, *pc) == NULL) {
		if (isspace(*pc)) {
			pc++;
			continue;
		}

		int token = emitexpr_lexer(pc, &pe);
		token_t new_tk = {0, token, 0, 0, 0, NULL};

		if (token == T_OP_L_PAREN || token == T_OP_L_SUBSTR) {
			int tk_prec = PRECEDENCE_OF(token);
			new_tk.precedence = tk_prec;
			stack_push(expr_stack, new_tk);

			ast_node_t *subexpr_root = emitexpr_parse(pc + 1, 
					(token == T_OP_L_PAREN) ? ")" : "]", 
					&pe, &redu_tk, pradix, pinit_radix,
					kmer_len, qmax);

			stack_push(expr_stack, redu_tk);
			if (subexpr_root) {
				stack_push(ast_stack, subexpr_root);
			}
		}
		else if (IS_OP(token)) {
			int tk_prec = PRECEDENCE_OF(token);
			if (tk_prec >= max_precedence) {
				max_precedence = tk_prec;
			}
			else {
				if (token != T_OP_R_SUBSTR && token != T_OP_R_PAREN) {
					max_precedence = emitexpr_reduce(&expr_stack, 
							&ast_stack, max_precedence, tk_prec, pradix, 
							pinit_radix, kmer_len, qmax);

					if (max_precedence < 0) {
						INFO("Failed to reduce expression!!\n");
						exit(1);
					}
				}
			}

			new_tk.precedence = tk_prec;
			stack_push(expr_stack, new_tk);

			if (token == T_OP_R_PAREN || token == T_OP_R_SUBSTR) {
				max_precedence = emitexpr_reduce(&expr_stack, 
						&ast_stack, PRECEDENCE_OF(token) + 1, 0, 
						pradix, pinit_radix, kmer_len, qmax);
			}
		}
		else {
			char *endptr;
			int num;
			int valid_func_flag;
			function_t const *f;

			switch (token) {
				case T_STR_LITERAL:
					valid_func_flag = 0;
					f = emitexpr_func_table;
					for (; f->name != NULL; f++) {
						if (strncmp(f->name, pc, pe-pc) == 0) {
							valid_func_flag = 1;
							break;
						}
					}

					if (valid_func_flag) {
						new_tk.token_type = T_OP_FUNC_NAME;
						new_tk.immediate = f->func_index;
						new_tk.token_as = f->return_type;
					}
					else {
						/* ERROR */
						INFO("Invalid keyword!\n");
						exit(1);
					}

					break;
				case T_NUM_LITERAL:
					endptr = NULL;
					num = strtol(pc, &endptr, 0);
					if (endptr != pe) {
						/* error */
						INFO("Invalid numeric constant!\n");
						exit(1);
					}

					new_tk.immediate = num;

					break;

				case T_KMER_LEN:
					new_tk.token_type = T_NUM_LITERAL;
					new_tk.immediate = kmer_len;

					break;

				case T_XT:
				case T_ST:
				case T_YT:
					/* equals to x[1..k] or s[1..k] */
					if (token == T_XT || token == T_ST)
						new_tk.token_as = T_BASES;

					new_tk.immediate = ((uint64_t) kmer_len << 32) | 1;

					break;
			}

			stack_push(expr_stack, new_tk);

			if (token != T_STR_LITERAL) {
				ast_node_t *_ast_node = malloc(sizeof(*_ast_node));
				_ast_node->lhs = NULL;
				_ast_node->rhs = NULL;
				_ast_node->token = new_tk;
				stack_push(ast_stack, _ast_node);
			}
		}

		pc = pe;
	}

	emitexpr_reduce(&expr_stack, &ast_stack, max_precedence, 0, 
			pradix, pinit_radix, kmer_len, qmax);

	/* pass out the position where (sub)-expression ends */
	if (endexpr) {
		*endexpr = pc;
	}

	/* pass out reduced token */
	if (expr_stack.used == 1) {
		if (reduced_token) {
			*reduced_token = expr_stack.data[0];
		}
		else {
			/* initial functional call instead of recursive call */
			token_t *tk = &expr_stack.data[0];
			if (tk->token_as == T_BASES) {
				int ctx_len = (tk->immediate >> 32) - 
							(tk->immediate & UINT32_MAX) + 1;

				if (tk->token_type == T_ST && 
						(tk->immediate >> 32) == kmer_len) {
					*pinit_radix = *pradix * (int) pow(4, ctx_len-1);
				}

				*pradix *= (int) pow(4, ctx_len);
			}
			else if (tk->token_type == T_COUNTER) {
				int ctx_len = (tk->immediate >> 32) - 
							(tk->immediate & UINT32_MAX) + 1;
				*pradix *= (ctx_len + 1);
			}
			else if (tk->token_type == T_INDICATOR) {
				*pradix *= 2;
			}
			else if (tk->token_type == T_YT) {
				int ctx_len = (tk->immediate >> 32) - 
							(tk->immediate & UINT32_MAX) + 1;

				*pradix *= (int) pow(qmax, ctx_len);
			}
			else {
				INFO("Final element can't be processed!.\n");
				exit(1);
			}

			INFO("[CTX] # of emission contexts: %zu\n", *pradix);
		}
	}
	else {
		INFO("(Sub)-expression can not be fully reduced.\n");
		exit(1);
	}

	/* flatten tree structure to an array */
	ast_node_t *rt_node = NULL;
	if (ast_stack.used) {
		rt_node = ast_stack.data[0];
	}

	stack_destroy(expr_stack);
	stack_destroy(ast_stack);

	return rt_node;
}

int emitexpr_lexer(char *pc, char **pe)
{
	char ch = *pc;
	*pe = pc + 1;

	switch (ch) {
		case '(':
			return T_OP_L_PAREN;
			break;
		case ')':
			return T_OP_R_PAREN;
			break;
		case '-':
			return T_OP_MINUS;
			break;
		case '+':
			return T_OP_PLUS;
			break;
		case '[':
			return T_OP_L_SUBSTR;
			break;
		case ']':
			return T_OP_R_SUBSTR;
			break;
		case '.':
			if (*(pc+1) == '.') {
				*pe = pc + 2;
				return T_OP_SLICE;
			}
			else {
				return T_INVALID;
			}
			break;
		case ',':
			return T_OP_COMMA;
			break;
		case '!':
			if (*(pc+1) == '=') {
				*pe = pc + 2;
				return T_OP_NE;
			}
			else {
				return T_INVALID;
			}
			break;
		case '<':
			if (*(pc+1) == '=') {
				*pe = pc + 2;
				return T_OP_LE;
			}
			else {
				return T_OP_LT;
			}
			break;
		case '>':
			if (*(pc+1) == '=') {
				*pe = pc + 2;
				return T_OP_GE;
			}
			else {
				return T_OP_GT;
			}
			break;
		case '=':
			return T_OP_EQ;
			break;
		default:
			if (isdigit(ch)) {
				char *p = pc + 1;
				for (; isdigit(*p); p++);
				*pe = p;
				return T_NUM_LITERAL;
			}
			else if (isalpha(ch)) {
				if (ch == 'x' && !isalpha(*(pc+1))) {
					return T_XT;
				}

				if (ch == 's' && !isalpha(*(pc+1))) {
					return T_ST;
				}

				if (ch == 'y' && !isalpha(*(pc+1))) {
					return T_YT;
				}

				if (ch == 'k' && !isalpha(*(pc+1))) {
					return T_KMER_LEN;
				}

				char *p = pc + 1;
				for (; isalpha(*p); p++);
				*pe = p;
				return T_STR_LITERAL;
			}
	}

	return T_INVALID;
}

int emitexpr_eval(oplist_t *oplst, uint64_t st, uint64_t xt, 
		const char *conv_q, int qmax)
{
	static uint64_t st_data[8];
	static stack_register_t st_int = {8, 0, st_data, UINT64_MAX};
	static int yt_len = 0;
	st_int.used = 0;

	uint64_t imm_reg_types[N_REGISTERS];
	uint64_t imm_registers[N_REGISTERS]; 
	uint64_t poped = 0;

	register int n_ops = oplst->n_ops;
	register op_t *op = oplst->ops;
	register int tmp_radix = 1;
	register int radix = 1;
	register int bctx = 1;

	for (int i = 0; i < n_ops; i++, op++) {
		switch (op->op_code) {
			case OP_LOAD_XT:
			case OP_LOAD_ST:
				imm_registers[op->op_dst] = op->op_code == OP_LOAD_ST ? st : xt;
				imm_reg_types[op->op_dst] = T_BASES;
				break;
			case OP_LOAD_YT:
				imm_registers[op->op_dst] = conv_q;
				imm_reg_types[op->op_dst] = T_YT;
				break;
			case OP_MOVE:
				imm_registers[op->op_dst] = (uint64_t) op->param;
				break;
			case OP_PUSH:
				stack_push(st_int, (uint64_t) op->param);
				break;
			case OP_POP:
				poped = stack_pop(st_int);
				imm_registers[op->op_dst] = poped;
				break;
			case OP_PEEK:
				poped = *stack_top(st_int);
				break;
			case OP_CMP_NE:
			case OP_CMP_EQ:
				do {
					uint64_t xor_sx = imm_registers[op->op_src1] ^
						imm_registers[op->op_src2];

					uint64_t xor_lo = xor_sx & KBITS_LO_MASK;
					uint64_t xor_hi = xor_sx & KBITS_HI_MASK;
					xor_lo |= (xor_hi >> 1);

					if (op->op_code == OP_CMP_EQ)
						xor_lo = ~xor_lo & KBITS_LO_MASK;

					imm_registers[op->op_dst] = xor_lo;
				} while(0);

				break;

			case OP_SUBSTRQ:
				do {
					uint32_t lb = ((uint64_t) op->param & UINT32_MAX) - 1;
					uint32_t len = ((uint64_t) op->param >> 32) - lb;
					uint64_t kid = imm_registers[op->op_src1];
					imm_registers[op->op_dst] = (kid >> (lb << 1)) & 
						((1UL << (len << 1)) - 1);
				} while (0);

				break;

			case OP_ADD:
				imm_registers[op->op_dst] = imm_registers[op->op_src1] + 
					(uint64_t) op->param;
				break;

			case OP_SET_YT_LEN:
				yt_len = (int) op->param;
				break;

			case OP_CMP_LT:
				op_cmp_tmpl(<);
				break;

			case OP_CMP_LE:
				op_cmp_tmpl(<=);
				break;

			case OP_CMP_GT:
				op_cmp_tmpl(>);
				break;

			case OP_CMP_GE:
				op_cmp_tmpl(>=);
				break;

			case OP_CALL:
				do {
					uint64_t fxn_idx = (uint64_t) op->param;
					uint64_t bools = imm_registers[op->op_src1];

					/* apply mask based on length of boolean vector */
					bools = bools & ((1 << (poped << 1)) - 1);
					if (fxn_idx == FXN_ALL) {
						bools = (_mm_popcnt_u64(bools) == poped);
					}
					else if (fxn_idx == FXN_ANY) {
						bools = bools > 0;
					}
					else if (fxn_idx == FXN_COUNT) {
						bools = _mm_popcnt_u64(bools);
					}

					imm_registers[op->op_dst] = bools;
				} while (0);

				break;
			case OP_TYPE:
				imm_reg_types[op->op_dst] = (uint64_t) op->param;
				break;

			case OP_RADIX_SET:
				tmp_radix = imm_registers[op->op_src1];
				break;

			case OP_REDUCE:
				do {
					uint64_t rtype = imm_reg_types[op->op_src2];
					if (rtype == T_BASES) {
						radix *= (int) pow(4, tmp_radix);
					}
					else if (rtype == T_INDICATOR) {
						radix *= 2;
					}
					else if (rtype == T_COUNTER) {
						radix *= (tmp_radix + 1);
					}
					else if (rtype == T_YT) {
						radix *= (int) pow(qmax, tmp_radix);
					}

					uint64_t imm1 = imm_reg_types[op->op_src1] == T_YT ?
						_eval_yt((char *) imm_registers[op->op_src1],
								yt_len ,qmax) : 
						imm_registers[op->op_src1];

					uint64_t imm2 = rtype == T_YT ?
						_eval_yt((char *) imm_registers[op->op_src2], 
								yt_len, qmax) : 
						imm_registers[op->op_src2];
					
					bctx = imm1 * radix + imm2;
					imm_registers[op->op_src1] = bctx;
				} while (0);

				break;
			case OP_OUTPUT:
				bctx = imm_registers[op->op_src1];
				break;
		}
	}

	return bctx;
}

void emitexpr_compile(ast_node_t *ast_root, oplist_t *oplst)
{
	emitexpr_ast_traverse(ast_root, 0, emitexpr_cb_compile, 
			(void *) oplst, 1);

	/* no reduction, directly output the context index */
	if (oplst->ops[oplst->n_ops-1].op_code != OP_REDUCE) {
		uint64_t imm_reg = stack_pop(*(oplst->sr));
		emitexpr_push_back_op(oplst, 
				OP_OUTPUT, 
				imm_reg, 0xFF, 0xFF, NULL);

		INFO("OUTPUT $%lu\n", imm_reg);
	}
}
