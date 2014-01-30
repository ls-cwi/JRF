#define __SEQUTIL_C__

#include "seqUtil.h"

typedef struct seqMem
{
	void *pos;
	struct seqMem *next;
} seqMem;

seqMem *start;
seqMem *current;

void seqMemInit()
{
	start = NULL;
}

void* seqMalloc(int size)
{
	if (start == NULL)
	{
		start = (seqMem*)malloc(sizeof(seqMem));
		memset(start, '\0', sizeof(seqMem));
		current = start;
	}
	else
	{
		current->next = (seqMem*)malloc(sizeof(seqMem));
		memset(current->next, '\0', sizeof(seqMem));
		current = current->next;
	}
	current->pos = malloc(size);
	memset(current->pos, '\0', size);
	return(current->pos);
}

void seqFreeAll()
{
	while (start != NULL)
	{
		current = start->next;
		free(start->pos);
		free(start);
		start = current;
	}

	start = NULL;
}

void seqFree(void* pos)
{
	seqMem *node, *prenode;
	node = start;
	prenode = start;
	while (node != NULL)
	{
		if (node->pos == pos)
		{
			free(node->pos);
			if (node == start)
			{
				start = node->next;
			}
			else if (node->next == NULL)
			{
				current = prenode;
				prenode->next = NULL;
			}
			else
			{
				prenode->next = node->next;
			}
			free(node);
			break;
		}

		prenode = node;
		node = node->next;
	}
}

void inputString(char *input, char **ppcStr, int *iLen, int *iMaxLen)
{
	int inputLen;
	char *temp;
	inputLen = strlen(input);
	if (inputLen == 0)
	{
		return;
	}
	while (*iMaxLen < (*iLen + inputLen) + 1)
	{
		*iMaxLen = *iMaxLen + APPEND_LEN;
	}
	temp = (char*)seqMalloc(*iMaxLen);
	if (*ppcStr == NULL)
	{
		memcpy(temp, input, inputLen);
	}
	else
	{
		memcpy(temp, *ppcStr, *iLen);
		strcat(temp, input);
	}
	*iLen = *iLen + inputLen;
	if (*ppcStr != NULL)
	{
		seqFree(*ppcStr);
	}
	*ppcStr = temp;
}

