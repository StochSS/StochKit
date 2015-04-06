#ifndef _btree_h_
#define _btree_h_

namespace STOCHKIT
{
	struct btree
	{
		std::size_t cellIndex;
		btree *leftChild;
		btree *rightChild;
		btree *parent;
		bool appendLeft;
		bool hasChild;
	};

	void addNode(btree *root, btree *node)
	{
		if(root!=node)
		{
			if(root->appendLeft)
			{
				if(root->leftChild==NULL)
				{
					root->leftChild=node;
					node->parent=root;
					root->hasChild=true;
				}
				else
					addNode(root->leftChild, node);
				root->appendLeft=false;
			}
			else
			{
				if(root->rightChild==NULL)
				{
					root->rightChild=node;
					node->parent=root;
					root->hasChild=true;
				}
				else
					addNode(root->rightChild, node);
				root->appendLeft=true;
			}
		}
	}

	void clearTree(btree tree[], std::size_t NumberOfSubvolumes)
	{
		for(std::size_t i=0; i<NumberOfSubvolumes; i++)
		{
			tree[i].cellIndex=NumberOfSubvolumes;
			tree[i].appendLeft=true;
		}
	}

	void initializeTree(btree *root, btree **nodeArray, double nextReactionTime[], std::size_t NumberOfSubvolumes)
	{
		std::size_t index, temp;
		btree *restore_root=root;
		nodeArray[0]=root;
		root->cellIndex=0;
		for(std::size_t i=1; i<NumberOfSubvolumes; i++)
		{
			root=restore_root;
			index=i;
			while(root->cellIndex!=NumberOfSubvolumes)
			{
				if(nextReactionTime[index]<nextReactionTime[root->cellIndex])
				{
					nodeArray[index]=root;
					temp=index;
					index=root->cellIndex;
					root->cellIndex=temp;
				}
				root->appendLeft=!root->appendLeft;//change appendLeft after adding the node (suppose to be done after we choose the child, we do it here since we will lost its pointer if we go to the child)
				root=root->appendLeft?root->rightChild:root->leftChild;//choose the child reverses the appendLeft since we changed appendLeft before we go to the child
			}
			//reach the node where we append the info
			nodeArray[index]=root;
			root->cellIndex=index;
		}
	}

	void maintain(btree *root, btree **nodeArray, double nextReactionTime[], bool downward=true)
	{
		std::size_t temp;
		btree *child;
		if(downward)
		{
			while(root->hasChild)
			{
				if(root->rightChild==NULL)
					child=root->leftChild;
				else if(root->leftChild==NULL)
					child=root->rightChild;
				else
					child=nextReactionTime[root->leftChild->cellIndex]<nextReactionTime[root->rightChild->cellIndex]?root->leftChild:root->rightChild;
				if(nextReactionTime[root->cellIndex]>nextReactionTime[child->cellIndex])
				{
					temp=root->cellIndex;
					root->cellIndex=child->cellIndex;
					child->cellIndex=temp;
					nodeArray[root->cellIndex]=root;
					nodeArray[child->cellIndex]=child;
					root=child;
				}
				else
					break;
			}
		}
		else
		{
			while(root->parent!=NULL)
			{
				if(nextReactionTime[root->cellIndex]<nextReactionTime[root->parent->cellIndex])
				{
					temp=root->cellIndex;
					root->cellIndex=root->parent->cellIndex;
					root->parent->cellIndex=temp;
					nodeArray[root->cellIndex]=root;
					nodeArray[root->parent->cellIndex]=root->parent;
					root=root->parent;
				}
				else
					break;
			}
		}
	}
}
#endif