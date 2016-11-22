#ifndef __PMR_EMITEXPR_H__
#define __PMR_EMITEXPR_H__

#include <stdlib.h>
#include <stdint.h>

#define N_REGISTERS 8

/*
 * Supported tokens:
 * operators and functions: -, _{ }, [ ], .., I( ), C( ), =, >, <, 
 * immediates: constant numerics, return value of I(.) and C(.)
 * nucleotides: x, y, s, and x[a..b] etc
 */

typedef struct ast_node_s ast_node_t;
typedef struct function_s function_t;
typedef struct stack_s stack_t;
typedef struct token_s token_t;
typedef struct op_s op_t;
typedef struct oplist_s oplist_t;

//extern const function_t emitexpr_func_table[];
extern const int emitexpr_precedence_table[];
extern const int emitexpr_assoc_table[];

#define IS_OP(t) ((t) < T_STR_LITERAL)
#define PRECEDENCE_OF(t) (emitexpr_precedence_table[(t)])
#define ASSOC_OF(t) (emitexpr_assoc_table[(t)])

#define op_cmp_tmpl(sym)									\
	do {													\
		uint64_t cmp_flag = 0;								\
		char *q = (char *) imm_registers[op->op_src1];		\
		uint64_t imm = imm_registers[op->op_src2];			\
		for (int i = 0; i < yt_len; ++i) {					\
			cmp_flag |= (q[i] sym imm) << (i << 1);			\
		}													\
		imm_registers[op->op_dst] = cmp_flag;				\
	} while(0);


ast_node_t *emitexpr_parse(char *expr, const char *end_chars, 
		char **endexpr, token_t *reduced_token, size_t *pradix, size_t *pinit_radix,
		int kmer_len, int qmax);
int emitexpr_eval(oplist_t *oplst, uint64_t st, uint64_t xt, 
		const char *conv_q, int qmax);
void emitexpr_compile(ast_node_t *ast_root, oplist_t *oplst);

/* FIXME: currently cannot handle cases with "N" bases */
//int emitexpr_eval(token_t *ast, kbits_t st, kbits_t xt);

enum {
	FXN_ALL, FXN_ANY, FXN_COUNT
};

enum {
	OP_RESET_RADIX,
	OP_LOAD_ST,
	OP_LOAD_XT,
	OP_LOAD_YT,
	OP_MOVE,
	/* substring on uint64_t */
	OP_SUBSTRQ,
	OP_ADD,
	OP_SET_YT_LEN,
	OP_CMP_EQ,
	OP_CMP_NE,
	OP_CMP_LT,
	OP_CMP_LE,
	OP_CMP_GT,
	OP_CMP_GE,
	OP_CALL,
	OP_PUSH,
	OP_POP,
	OP_PEEK,
	OP_REDUCE,
	OP_RADIX_SET,
	OP_OUTPUT,
	OP_TYPE,
	OP_INVALID = 0xFF
};

enum {
	/* [, ] */
	T_OP_L_SUBSTR, T_OP_R_SUBSTR,
	/* func(, (, ) */
	T_OP_FUNC_NAME, T_OP_L_PAREN, T_OP_R_PAREN,
	T_OP_EQ, T_OP_LT, T_OP_LE, T_OP_GT, T_OP_GE, 
	/* != */
	T_OP_NE,
	/* +, - */
	T_OP_PLUS, T_OP_MINUS,
	/* .. */
	T_OP_SLICE, 
	/* , */
	T_OP_COMMA,
	T_OP_SUBSTR_OF, /* [a..b], [k-2..k] */
	T_OP_FUNCTION, /* all, any, count */

	T_STR_LITERAL,
	T_NUM_LITERAL, /* 1, 2, k ... */
	T_SLICE, /* a..b, k-2..k */
	T_XT, /* x, x[a..b] */
	T_ST, /* s, s[a..b] */
	T_YT, /* y, y[a..b] */
	T_BASES, /* T_XT, T_ST */
	T_KMER_LEN, /* k */
	T_BOOLS, 
	T_INDICATOR, /* any(a=b), all(x[k-3..k-1] != s[k-3..k-1]) ... */
	T_COUNTER,  /* count(a=b) */
	T_INVALID = 0xFF
};

enum {
	ASSOC_NONE, ASSOC_BINARY, ASSOC_UNITARY_RIGHT, ASSOC_UNITARY_LEFT
};

struct function_s {
	const char *name;
	int func_index;
	int arg_type;
	int return_type;
};

struct token_s {
	uint32_t precedence : 4;
	uint32_t token_type : 8;
	/* for a function token, it can be treated "as" if a token of which
	 * the type is its returned type. */
	uint32_t token_as : 8;
	uint32_t reserved : 12;
	/* if a token has a immediate value that can be fit into a 64-bit
	 * type, it can be put into `immediate` */
	uint64_t immediate;
	void *data;
};

struct ast_node_s {
	token_t token;
	ast_node_t *lhs;
	ast_node_t *rhs;
};

#define stack_typedef(type) \
	typedef struct { \
		int size;			\
		int used;			\
		type *data;			\
		type dummy;			\
	}

#define stack_init(type, st, sz, dm)			\
	do {										\
		(st).size = sz;							\
		(st).used = 0;							\
		(st).data = malloc(sz * sizeof(type));	\
		(st).dummy = dm;						\
	} while(0);

#define stack_destroy(s)						\
	do {										\
		free((s).data);							\
	} while(0);									\

#define stack_empty(s) ((s).used == 0)
#define stack_top(s) (&((s).data[(s).used-1]))
#define stack_bottom(s) (&((s).data[0]))

#define stack_push(s, t)				\
	do {								\
		if ((s).used < (s).size) {		\
			(s).data[(s).used++] = t;	\
		}								\
	} while (0);

#define stack_pop(s) \
	((s).used == 0) ? (s).dummy : (s).data[--(s).used]

stack_typedef(token_t) stack_token_t;
stack_typedef(ast_node_t *) stack_ast_t;
stack_typedef(uint64_t) stack_register_t;

struct op_s {
	uint32_t op_code : 8;
	uint32_t op_src1 : 8;
	uint32_t op_src2 : 8;
	uint32_t op_dst : 8;
	void *param;
};

struct oplist_s {
	stack_register_t *sr;
	int n_ops;
	op_t *ops;
};


#endif
