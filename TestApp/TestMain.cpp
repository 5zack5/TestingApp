#include <iostream>
#include <ostream>
#include <vector>
#include <deque>

using namespace std;

int TOH(int source, int destination, int intermediate, int numDisks);
int ContigousSum(int *array, int numElements);
int FindMissing(int *elem, int numElements);
int SubArraySum(int *elem, int numElements, int targetSum);
int Sort012(int *elem, int numElements);
int FindEq(int *elem, int numElements);
int MaxSumSubSequence(int *elem, int numElements);
int LeadersInArray(int *elem, int numElements);
int MaxProfit(int *elem, int numElements);
int MaxSubArrayK(int *elem, int numElements, int k);

int main(int argc, char *argv[])
{
	//TOH
	//int numDisks     = atoi(argv[1]);
	//TOH(1,3,2,numDisks);

	//ContigousSum
	//int numElements = atoi(argv[1]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//   elements[i] = atoi(argv[1+i+1]);
	//ContigousSum(elements,numElements);

	//Missing number
	//int numElements = atoi(argv[1]);
	//int *elements = new int[numElements-1];
	//for (int i=0;i<numElements-1;i++)
	//   elements[i] = atoi(argv[1+i+1]);
	//FindMissing(elements,numElements);
	
	//Subarray with specified sum
	//int numElements = atoi(argv[1]);
	//int targetSum = atoi(argv[2]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//   elements[i] = atoi(argv[3+i]);
	//SubArraySum(elements,numElements, targetSum);
	

	////Sort 012 array
	//int numElements = atoi(argv[1]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//   elements[i] = atoi(argv[1+i+1]);
	//Sort012(elements,numElements);
	

	//Equilibrium point of an array
	//int numElements = atoi(argv[1]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//   elements[i] = atoi(argv[1+i+1]);
	//FindEq(elements,numElements);

	//delete [] elements;


	//MaxSumSubsequence
	//int numElements = atoi(argv[1]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//   elements[i] = atoi(argv[1+i+1]);
	//MaxSumSubSequence(elements,numElements);

	//delete [] elements;
	
	//Leaders In Array
	//int numElements = atoi(argv[1]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//   elements[i] = atoi(argv[1+i+1]);
	//LeadersInArray(elements,numElements);

	//Max Profit
	//int numElements = atoi(argv[1]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//   elements[i] = atoi(argv[1+i+1]);
	//MaxProfit(elements,numElements);

	//delete [] elements;


	////Max subarray K
	int numElements = atoi(argv[1]);
	int *elements = new int[numElements];
	for (int i=0;i<numElements;i++)
	   elements[i] = atoi(argv[1+i+1]);
	MaxSubArrayK(elements,numElements);

	//delete [] elements;
}

int MaxProfit(int *elem, int numElements)
{
	int globalMin = elem[0];
	int MaxGain = -100000;

	for (int i=1;i<numElements;i++)
	{
		if ((elem[i]-globalMin) > MaxGain)
			MaxGain = elem[i]-globalMin;

		if (elem[i]<globalMin)
			globalMin = elem[i];

	}
	cout << "Max gain:" << MaxGain << endl;
	return 0;
}

int MaxSubArrayK(int *elem, int numElements, int k)
{
	deque<int> deQue(k);


	return 0;
}

int LeadersInArray(int *elem, int numElements)
{
  std::vector<int> indices;
  indices.push_back(numElements-1);

  int maxSoFar=elem[numElements-1];
  for (int i=numElements-2;i>=0;i--)
  {
	  maxSoFar = max(maxSoFar,elem[i]);
	  if (maxSoFar==elem[i])
		  indices.push_back(i);
  }

  cout << "Leaders in the array:";
  for(int i=indices.size()-1;i>=0;i--)
	  cout  << elem[indices[i]] << " " ;
  cout << endl;

	return 0;
}

int MaxSumSubSequence(int *elem, int numElements)
{
	int globalMaxSum=elem[0];
	int sum=elem[0];
	int prev= elem[0];
	for (int i=1;i<numElements;i++)
	{
		if(elem[i]>prev)
			sum+=elem[i];
		else
			sum=elem[i];

		if (sum>globalMaxSum)
			globalMaxSum = sum;

		prev = elem[i];
	}

	cout << "Global Max Sum: " << globalMaxSum << endl;



	return 0;
}

int FindEq(int *elem, int numElements)
{
	int sum=0;
	for (int i=0;i<numElements;i++)
		sum+=elem[i];

	cout << sum << endl;
	int leftSum=0;
	for (int i=0;i<numElements-1;i++)
	{
		sum-=elem[i];
		cout << leftSum << ":" << sum-elem[i] << endl;
		if (leftSum==sum)
		{
			cout << "Equlibrium at index:" << i;
			break;
		}
		leftSum+=elem[i];
	}

	return 0;
}

void swap(int *elem, int index1, int index2)
{
   int temp = elem[index1];
   elem[index1] = elem[index2];
   elem[index2] = temp;
}

int Sort012(int *elem, int numElements)
{
	int lo=0;
	int md=0;
	int hi=numElements-1;

	while (md<=hi)
	{
		switch (elem[md])
		{
			case 0:
				swap(elem,lo++,md++);
				break;
			case 1:
				md++;
				break;
			case 2:
				swap(elem,md,hi--);
				break;
		}
	}

	cout << "Sorted Array:"<< endl;
	for (int i=0;i<numElements;i++)
		cout << elem[i] << " " ;
	cout << endl;

	return 0;

}

int SubArraySum(int *elem, int numElements, int targetSum)
{
	if (numElements==0)
	{
		cout << "No elements in the array!" << endl;
		return -1;
	}

	bool sumFound = false;
	int startIndex = 0;
	int stopIndex  = 0;
	int globalSum=elem[0];
	int sum=elem[0];
	

	for (int i=startIndex+1;i<numElements;i++)
	{
		sum+=elem[i];
		stopIndex = i;
		
		if (sum == targetSum)
		{
			sumFound  = true;
			stopIndex = i;
			break;
		}
		
		while (sum > targetSum)
		{
			sum-=elem[startIndex++];
			if (sum==targetSum)
			{
				sumFound = true;
				break;
			}
		}

	}
	
	if (sumFound)
		cout << "Subarray found! " << startIndex+1 << "," << stopIndex+1 << endl;
	else
		cout << "Subarray not found!" << endl;

	return 0;
}


int FindMissing(int *elem, int numElements)
{
	if (numElements==0)
	{
		cout << "No elements in the array!" << endl;
		return -1;
	}

	int expSum= numElements*(numElements+1)/2;
	int sum =0;
	for (int i=0;i<numElements-1;i++)
		sum+=elem[i];

	cout << "Missing number:" << expSum-sum << endl;
	return 0;

}

int ContigousSum(int *elem, int numElements)
{
	if (numElements==0)
	{
		cout << "No elements in the array!" << endl;
		return -1;
	}

	int maxSum=elem[0];
	int sum=elem[0];

	for (int i=1;i<numElements;i++)
	{
		sum = max(elem[i],sum+elem[i]);
		if (sum > maxSum)
			maxSum = sum;
	}

	cout << "Max contiguous sum :" << maxSum;
	return 0;
}

int TOH(int source, int destination, int intermediate, int numDisks)
{
	if (numDisks==0)
		return 0;

	TOH(source, intermediate, destination, numDisks-1);
	cout << "Moving disk "  << numDisks << " from rod " << source << " to rod " << destination << endl;
	TOH(intermediate,destination,source, numDisks-1);
	

	return 0;
}



