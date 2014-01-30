#ifndef __NEWICKFORM_H__
#define __NEWICKFORM_H__

typedef struct newick_child
{
	struct newick_node *node;
	struct newick_child *next;
} newick_child;

typedef struct newick_node
{
	char *taxon;
	char *seq;
	float dist;
	int childNum;
	struct newick_child *child;
	struct newick_node *parent;
} newick_node;

//short sDEBUG;

//#ifdef __NEWICKFORM_C__
newick_node* parseTree(char *str);
void printTree(newick_node *root);

//#else
//extern newick_node* parseTree(char *str);
//extern void printTree(newick_node *root);
//
//#endif

#endif

