/**
 *  @file Huffman.c
 *  @author Sheng Di
 *  @date Aug., 2016
 *  @brief Customized Huffman Encoding, Compression and Decompression functions
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Huffman.h"
#include "sz.h"

int stateNum;
int allNodes;

//struct node_t pool[allNodes] = {{0}};
node pool;
node *qqq = NULL;
node *qq = NULL;
int n_nodes = 0;
int qend = 1;
unsigned long **code = NULL;//TODO
unsigned char *cout = NULL;
int n_inode = 0;
 
node new_node(int freq, unsigned int c, node a, node b)
{
	node n = pool + n_nodes++;
	if (freq) 
	{
		n->c = c;
		n->freq = freq;
		n->t = 1;
	}
	else {
		n->left = a; 
		n->right = b;
		n->freq = a->freq + b->freq;
		n->t = 0;
		//n->c = 0;
	}
	return n;
}
 
node new_node2(unsigned int c, unsigned char t)
{
	pool[n_nodes].c = c;
	pool[n_nodes].t = t;
	return pool+n_nodes++;
} 
 
/* priority queue */
void qinsert(node n)
{
	int j, i = qend++;
	while ((j = i / 2)) {
		if (qq[j]->freq <= n->freq) break;
		qq[i] = qq[j], i = j;
	}
	qq[i] = n;
}
 
node qremove()
{
	int i, l;
	node n = qq[i = 1];
 
	if (qend < 2) return 0;
	qend--;
	while ((l = i * 2) < qend) {
		if (l + 1 < qend && qq[l + 1]->freq < qq[l]->freq) l++;
		qq[i] = qq[l], i = l;
	}
	qq[i] = qq[qend];
	return n;
}
 
/* walk the tree and put 0s and 1s */
/**
 * @out1 should be set to 0.
 * @out2 should be 0 as well.
 * @index: the index of the byte
 * */
void build_code(node n, int len, unsigned long out1, unsigned long out2)
{
	if (n->t) {
		code[n->c] = (unsigned long*)malloc(2*sizeof(unsigned long));
		if(len<=64)
		{
			(code[n->c])[0] = out1 << (64 - len);
			(code[n->c])[1] = out2;
		}
		else
		{
			(code[n->c])[0] = out1;
			(code[n->c])[1] = out2 << (128 - len);
		}
		cout[n->c] = (unsigned char)len;
		return;
	}
	int index = len >> 6; //=len/64
	if(index == 0)
	{
		out1 = out1 << 1;
		out1 = out1 | 0;
		build_code(n->left, len + 1, out1, 0);
		out1 = out1 | 1;
		build_code(n->right, len + 1, out1, 0);		
	}
	else
	{
		if(len%64!=0)
			out2 = out2 << 1;
		out2 = out2 | 0;
		build_code(n->left, len + 1, out1, out2);
		out2 = out2 | 1;
		build_code(n->right, len + 1, out1, out2);	
	}
}

void init(int *s, int length)
{
	int i, *freq = (int *)malloc(allNodes*sizeof(int));
	int index;
	memset(freq, 0, allNodes*sizeof(int));
	for(i = 0;i < length;i++) 
	{
		//index = 0;
		//index = (index | s[i])<<8;
		//index = index | s[i+1];
		index = s[i];
		freq[index]++;		
	}
 
	for (i = 0; i < allNodes; i++)
		if (freq[i]) 
			qinsert(new_node(freq[i], i, 0, 0));
 
	while (qend > 2) 
		qinsert(new_node(0, 0, qremove(), qremove()));
 
	build_code(qq[1], 0, 0, 0);
	free(freq);
}
 
void encode(int *s, int length, unsigned char *out, int *outSize)
{
	int i = 0;
	unsigned char curByte, bitSize = 0, byteSize, byteSizep;
	int state;
	unsigned char *p = out;
	int lackBits = 0;
	//long totalBitSize = 0, maxBitSize = 0, bitSize21 = 0, bitSize32 = 0;
	for (i = 0;i<length;i++) 
	{
		//state = 0;
		//state = (state | s[i])<<8;
		//state = state | s[i+1];
		
		state = s[i];
		bitSize = cout[state];	
		
		//printf("%d %d : %d %u\n",i, state, bitSize, (code[state])[0] >> (64-cout[state])); 
		//debug: compute the average bitSize and the count that is over 32... 	
		/*if(bitSize>=21)
			bitSize21++;
		if(bitSize>=32)
			bitSize32++;
		if(maxBitSize<bitSize)
			maxBitSize = bitSize;
		totalBitSize+=bitSize;*/

		if(lackBits==0)
		{
			byteSize = bitSize%8==0 ? bitSize/8 : bitSize/8+1; //it's equal to the number of bytes involved (for *outSize)
			byteSizep = bitSize/8; //it's used to move the pointer p for next data
			if(byteSize<=8)				
			{
				longToBytes_bigEndian(p, (code[state])[0]);
				p += byteSizep;
			}
			else //byteSize>8
			{
				longToBytes_bigEndian(p, (code[state])[0]);
				p += 8;			
				longToBytes_bigEndian(p, (code[state])[1]);
				p += (byteSizep - 8);		
			}
			*outSize += byteSize;
			lackBits = bitSize%8==0 ? 0 : 8 - bitSize%8;
		}
		else
		{
			*p = (*p) | (unsigned char)((code[state])[0] >> (64 - lackBits));			
			if(lackBits < bitSize)
			{
				p++;
				//(*outSize)++;
				long newCode = (code[state])[0] << lackBits;
				longToBytes_bigEndian(p, newCode);				

				if(bitSize<=64)
				{
					bitSize -= lackBits;
					byteSize = bitSize%8==0 ? bitSize/8 : bitSize/8+1;
					byteSizep = bitSize/8;
					p += byteSizep;
					(*outSize)+=byteSize;
					lackBits = bitSize%8==0 ? 0 : 8 - bitSize%8;
				}
				else //bitSize > 64
				{
					byteSizep = 7; //must be 7 bytes, because lackBits!=0
					p+=byteSizep;
					(*outSize)+=byteSize;
					
					bitSize -= 64;
					if(lackBits < bitSize)
					{
						*p = (*p) | (unsigned char)((code[state])[0] >> (64 - lackBits));
						p++;
						//(*outSize)++;						
						newCode = (code[state])[1] << lackBits;
						longToBytes_bigEndian(p, newCode);
						bitSize -= lackBits;
						byteSize = bitSize%8==0 ? bitSize/8 : bitSize/8+1;
						byteSizep = bitSize/8;
						p += byteSizep;
						(*outSize)+=byteSize;
						lackBits = bitSize%8==0 ? 0 : 8 - bitSize%8;						
					}
					else //lackBits >= bitSize
					{
						*p = (*p) | (unsigned char)((code[state])[0] >> (64 - bitSize));
						lackBits -= bitSize;
					}		
				}
			}
			else //lackBits >= bitSize
			{
				lackBits -= bitSize;
				if(lackBits==0)
					p++;
			}
		}
	}
//	for(i=0;i<stateNum;i++)
//		if(code[i]!=NULL) free(code[i]);
	/*printf("max bitsize = %d\n", maxBitSize);
	printf("bitSize21 ratio = %f\n", ((float)bitSize21)/length);
	printf("bitSize32 ratio = %f\n", ((float)bitSize32)/length);
	printf("avg bit size = %f\n", ((float)totalBitSize)/length);*/
}
 
void decode(unsigned char *s, int targetLength, node t, int *out)
{
	unsigned long i = 0, byteIndex = 0;
	int r, count=0;
	node n = t;
	char byte;
	for(i=0;count<targetLength;i++)
	{
		//if(count==2035)
		//	printf("%d\n",i);
		
		byteIndex = i>>3; //i/8
		r = i%8;
		if(((s[byteIndex] >> (7-r)) & 0x01) == 0)
			n = n->left;
		else
			n = n->right;
		if (n->t) {
			//putchar(n->c); 
			out[count] = n->c;
			n = t; 
			count++;
		}
	}
//	putchar('\n');
	if (t != n) printf("garbage input\n");
	return;
} 
	 
void pad_tree_uchar(unsigned char* L, unsigned char* R, unsigned int* C, unsigned char* t, unsigned int i, node root)
{
	C[i] = root->c;
	t[i] = root->t;
	node lroot = root->left;
	if(lroot!=0)
	{
		n_inode++;
		L[i] = n_inode;
		pad_tree_uchar(L,R,C,t, n_inode, lroot);
	}
	node rroot = root->right;
	if(rroot!=0)
	{
		n_inode++;
		R[i] = n_inode;
		pad_tree_uchar(L,R,C,t, n_inode, rroot);
	}
}  

void pad_tree_ushort(unsigned short* L, unsigned short* R, unsigned int* C, unsigned char* t, unsigned int i, node root)
{
	C[i] = root->c;
	t[i] = root->t;
	node lroot = root->left;
	if(lroot!=0)
	{
		n_inode++;
		L[i] = n_inode;
		pad_tree_ushort(L,R,C,t,n_inode, lroot);
	}
	node rroot = root->right;
	if(rroot!=0)
	{
		n_inode++;
		R[i] = n_inode;
		pad_tree_ushort(L,R,C,t,n_inode, rroot);
	}	
}

void pad_tree_uint(unsigned int* L, unsigned int* R, int* C, unsigned char* t, unsigned int i, node root)
{
	C[i] = root->c;
	t[i] = root->t;
	node lroot = root->left;
	if(lroot!=0)
	{
		n_inode++;
		L[i] = n_inode;
		pad_tree_uint(L,R,C,t,n_inode, lroot);
	}
	node rroot = root->right;
	if(rroot!=0)
	{
		n_inode++;
		R[i] = n_inode;
		pad_tree_uint(L,R,C,t,n_inode, rroot);
	}
}
 
unsigned int convert_HuffTree_to_bytes_anyStates(int nodeCount, unsigned char** out) 
{
	//printf("nodeCount=%d\n", nodeCount);
	if(nodeCount<=256)
	{
		unsigned char* L = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
		memset(L, 0, nodeCount*sizeof(unsigned char));
		unsigned char* R = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
		memset(R, 0, nodeCount*sizeof(unsigned char));
		unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
		memset(C, 0, nodeCount*sizeof(unsigned int));
		unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
		memset(t, 0, nodeCount*sizeof(unsigned char));

		pad_tree_uchar(L,R,C,t,0,qq[1]);
		
		unsigned int totalSize = 1+3*nodeCount*sizeof(unsigned char)+nodeCount*sizeof(unsigned int);	
		*out = (unsigned char*)malloc(totalSize*sizeof(unsigned char));
		(*out)[0] = (unsigned char)sysEndianType;
		memcpy(*out+1, L, nodeCount*sizeof(unsigned char));
		memcpy((*out)+1+nodeCount*sizeof(unsigned char),R,nodeCount*sizeof(unsigned char));
		memcpy((*out)+1+2*nodeCount*sizeof(unsigned char),C,nodeCount*sizeof(unsigned int));
		memcpy((*out)+1+2*nodeCount*sizeof(unsigned char)+nodeCount*sizeof(unsigned int), t, nodeCount*sizeof(unsigned char));
		free(L);
		free(R);
		free(C);
		free(t);
		return totalSize;

	}
	else if(nodeCount<=65536)
	{
		unsigned short* L = (unsigned short*)malloc(nodeCount*sizeof(unsigned short));
		memset(L, 0, nodeCount*sizeof(unsigned short));
		unsigned short* R = (unsigned short*)malloc(nodeCount*sizeof(unsigned short));
		memset(R, 0, nodeCount*sizeof(unsigned short));
		unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));	
		memset(C, 0, nodeCount*sizeof(unsigned int));		
		unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
		memset(t, 0, nodeCount*sizeof(unsigned char));		
		pad_tree_ushort(L,R,C,t,0,qq[1]);
		unsigned int totalSize = 1+2*nodeCount*sizeof(unsigned short)+nodeCount*sizeof(unsigned char) + nodeCount*sizeof(unsigned int);
		*out = (char*)malloc(totalSize);
		(*out)[0] = (unsigned char)sysEndianType;		
		memcpy(*out+1, L, nodeCount*sizeof(unsigned short));
		memcpy((*out)+1+nodeCount*sizeof(unsigned short),R,nodeCount*sizeof(unsigned short));
		memcpy((*out)+1+2*nodeCount*sizeof(unsigned short),C,nodeCount*sizeof(unsigned int));
		memcpy((*out)+1+2*nodeCount*sizeof(unsigned short)+nodeCount*sizeof(unsigned int),t,nodeCount*sizeof(unsigned char));
		free(L);
		free(R);
		free(C);
		free(t);		
		return totalSize;
	}
	else //nodeCount>65536
	{
		unsigned int* L = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
		memset(L, 0, nodeCount*sizeof(unsigned int));
		unsigned int* R = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
		memset(R, 0, nodeCount*sizeof(unsigned int));
		unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));	
		memset(C, 0, nodeCount*sizeof(unsigned int));
		unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
		memset(t, 0, nodeCount*sizeof(unsigned char));
		pad_tree_uint(L,R,C,t,0,qq[1]);
		unsigned int totalSize = 1+3*nodeCount*sizeof(unsigned int)+nodeCount*sizeof(unsigned char);
		*out = (unsigned char*)malloc(totalSize);
		(*out)[0] = (unsigned char)sysEndianType;
		memcpy(*out+1, L, nodeCount*sizeof(unsigned int));
		memcpy((*out)+1+nodeCount*sizeof(unsigned int),R,nodeCount*sizeof(unsigned int));
		memcpy((*out)+1+2*nodeCount*sizeof(unsigned int),C,nodeCount*sizeof(unsigned int));
		memcpy((*out)+1+3*nodeCount*sizeof(unsigned int),t,nodeCount*sizeof(unsigned char));
		free(L);
		free(R);
		free(C);
		free(t);
		return totalSize;		
	}
}

void unpad_tree_uchar(unsigned char* L, unsigned char* R, unsigned int* C, unsigned char *t, unsigned int i, node root)
{
	//root->c = C[i];
	if(root->t==0)
	{
		unsigned char l, r;
		l = L[i];
		if(l!=0)
		{
			node lroot = new_node2(C[l],t[l]);
			root->left = lroot;
			unpad_tree_uchar(L,R,C,t,l,lroot);
		}
		r = R[i];
		if(r!=0)
		{
			node rroot = new_node2(C[r],t[r]);
			root->right = rroot;
			unpad_tree_uchar(L,R,C,t,r,rroot);
		}
	}
}

void unpad_tree_ushort(unsigned short* L, unsigned short* R, unsigned int* C, unsigned char* t, unsigned int i, node root)
{
	//root->c = C[i];
	if(root->t==0)
	{
		unsigned short l, r;
		l = L[i];
		if(l!=0)
		{
			node lroot = new_node2(C[l],t[l]);
			root->left = lroot;
			unpad_tree_ushort(L,R,C,t,l,lroot);
		}
		r = R[i];
		if(r!=0)
		{
			node rroot = new_node2(C[r],t[r]);
			root->right = rroot;
			unpad_tree_ushort(L,R,C,t,r,rroot);
		}
	}
}

void unpad_tree_uint(unsigned int* L, unsigned int* R, int* C, unsigned char* t, unsigned int i, node root)
{
	//root->c = C[i];
	if(root->t==0)
	{
		unsigned int l, r;
		l = L[i];
		if(l!=0)
		{
			node lroot = new_node2(C[l],t[l]);
			root->left = lroot;
			unpad_tree_uint(L,R,C,t,l,lroot);
		}
		r = R[i];
		if(r!=0)
		{
			node rroot = new_node2(C[r],t[r]);
			root->right = rroot;
			unpad_tree_uint(L,R,C,t,r,rroot);
		}
	}
}

node reconstruct_HuffTree_from_bytes_anyStates(char* bytes, int nodeCount)
{
	//printf("nodeCount=%d\n", nodeCount);
	if(nodeCount<=256)
	{
		unsigned char* L = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
		unsigned char* R = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
		unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
		unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));	
		unsigned char cmpSysEndianType = bytes[0];
		if(cmpSysEndianType!=(unsigned char)sysEndianType)
		{
			unsigned char* p = bytes+1+2*nodeCount*sizeof(unsigned char);
			int i = 0, size = nodeCount*sizeof(unsigned int);
			while(1)
			{
				symTransform_4bytes(p);
				i+=sizeof(unsigned int);
				if(i<size)
					p+=sizeof(unsigned int);
				else
					break;
			}		
		}
		memcpy(L, bytes+1, nodeCount*sizeof(unsigned char));
		memcpy(R, bytes+1+nodeCount*sizeof(unsigned char), nodeCount*sizeof(unsigned char));
		memcpy(C, bytes+1+2*nodeCount*sizeof(unsigned char), nodeCount*sizeof(unsigned int));	
		memcpy(t, bytes+1+2*nodeCount*sizeof(unsigned char)+nodeCount*sizeof(unsigned int), nodeCount*sizeof(unsigned char));
		node root = new_node2(0,0);
		unpad_tree_uchar(L,R,C,t,0,root);
		free(L);
		free(R);
		free(C);
		free(t);
		return root;
	}
	else if(nodeCount<=65536)
	{
		unsigned short* L = (unsigned short*)malloc(nodeCount*sizeof(unsigned short));
		unsigned short* R = (unsigned short*)malloc(nodeCount*sizeof(unsigned short));
		unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));	
		unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
		unsigned char cmpSysEndianType = bytes[0];	
		if(cmpSysEndianType!=(unsigned char)sysEndianType)
		{
			unsigned char* p = bytes+1;
			int i = 0, size = 3*nodeCount*sizeof(unsigned int);
			while(1)
			{
				symTransform_4bytes(p);
				i+=sizeof(unsigned int);
				if(i<size)
					p+=sizeof(unsigned int);
				else
					break;
			}		
		}

		memcpy(L, bytes+1, nodeCount*sizeof(unsigned short));
		memcpy(R, bytes+1+nodeCount*sizeof(unsigned short), nodeCount*sizeof(unsigned short));
		memcpy(C, bytes+1+2*nodeCount*sizeof(unsigned short), nodeCount*sizeof(unsigned int));	

		memcpy(t, bytes+1+2*nodeCount*sizeof(unsigned short)+nodeCount*sizeof(unsigned int), nodeCount*sizeof(unsigned char));	

		node root = new_node2(0,0);
		unpad_tree_ushort(L,R,C,t,0,root);
		free(L);
		free(R);
		free(C);
		free(t);		
		return root;				
	}
	else //nodeCount>65536
	{
		unsigned int* L = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
		unsigned int* R = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
		unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));	
		unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));	
		unsigned char cmpSysEndianType = bytes[0];
		if(cmpSysEndianType!=(unsigned char)sysEndianType)
		{
			unsigned char* p = bytes+1;
			int i = 0, size = 3*nodeCount*sizeof(unsigned int);
			while(1)
			{
				symTransform_4bytes(p);
				i+=sizeof(unsigned int);
				if(i<size)
					p+=sizeof(unsigned int);
				else
					break;
			}
		}

		memcpy(L, bytes+1, nodeCount*sizeof(unsigned int));
		memcpy(R, bytes+1+nodeCount*sizeof(unsigned int), nodeCount*sizeof(unsigned int));
		memcpy(C, bytes+1+2*nodeCount*sizeof(unsigned int), nodeCount*sizeof(unsigned int));	
	
		memcpy(t, bytes+1+3*nodeCount*sizeof(unsigned int), nodeCount*sizeof(unsigned char));			
					
		node root = new_node2(0,0);
		unpad_tree_uint(L,R,C,t,0,root);
		free(L);
		free(R);
		free(C);
		free(t);
		return root;
	}
}

void encode_withTree(int *s, int length, unsigned char **out, int *outSize)
{
	int i, nodeCount = 0;
	unsigned char *treeBytes, buffer[4];
	
	init(s, length);
	for (i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);
	//printf("treeByteSize=%d\n", treeByteSize);
	*out = (unsigned char*)malloc(length*sizeof(int)+treeByteSize);
	intToBytes_bigEndian(buffer, nodeCount);
	memcpy(*out, buffer, 4);
	memcpy(*out+4, treeBytes, treeByteSize);
	free(treeBytes);
	int enCodeSize = 0;
	encode(s, length, *out+4+treeByteSize, &enCodeSize);
	*outSize = 4+treeByteSize+enCodeSize;
	
	//unsigned short state[length];
	//decode(*out+4+treeByteSize, enCodeSize, qqq[0], state);
	//printf("dataSeriesLength=%d",length );
}

/**
 * @par *out rememmber to allocate targetLength short_type data for it beforehand.
 * 
 * */
void decode_withTree(unsigned char *s, int targetLength, int *out)
{
	int encodeStartIndex;
	int nodeCount = bytesToInt_bigEndian(s);
	node root = reconstruct_HuffTree_from_bytes_anyStates(s+4, nodeCount);
	
	//sdi: Debug
/*	build_code(root, 0, 0, 0);
	int i;
	unsigned long code_1, code_2;
	for (i = 0; i < stateNum; i++)
		if (code[i])
		{		
			printf("%d: %lu,%lu ; %u\n", i, (code[i])[0],(code[i])[1], cout[i]);
			//code_1 = (code[i])[0];
		}*/
	
	if(nodeCount<=256)
		encodeStartIndex = 1+3*nodeCount*sizeof(unsigned char)+nodeCount*sizeof(unsigned int);
	else if(nodeCount<=65536)
		encodeStartIndex = 1+2*nodeCount*sizeof(unsigned short)+nodeCount*sizeof(unsigned char)+nodeCount*sizeof(unsigned int);
	else
		encodeStartIndex = 1+3*nodeCount*sizeof(unsigned int)+nodeCount*sizeof(unsigned char);
	decode(s+4+encodeStartIndex, targetLength, root, out);
}
