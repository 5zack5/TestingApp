#include <iostream>
#include <ostream>
#include <vector>
#include <deque>
#include <math.h>
#include <algorithm>
#include <queue>
#include <random>

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
//int MaxSubArrayK(int *elem, int numElements, int k);
int TaskOptimization(int *elem, int numTasks, int maxTime);
int MultiplyWithoutMult(int n1, int n2);
int KthSmallestElement(int *elem, int numElements, int k);
int RainWaterTrap(int *elem, int numElements);
int CoinSum(int *elem, int numElements, int targetSum);
int MedianLinearTime(int *elem, int numElements);
int MultiplyWithoutMultAdd(unsigned int &n1, unsigned int &n2);

int main(int argc, char *argv[])
{

	if (argc < 2)
		return 0;

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

	///Max subarray K
	//int numElements = atoi(argv[1]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//   elements[i] = atoi(argv[1+i+1]);
	//MaxSubArrayK(elements,numElements);

	//TaskOptimization
	//int numTasks = atoi(argv[1]);
	//int maxTime  = atoi(argv[2]);
	//int *elements = new int[numTasks];
	//for (int i=0;i<numTasks;i++)
	//   elements[i] = atoi(argv[3+i]);
	//TaskOptimization(elements,numTasks,maxTime);

	//int n1 = atoi(argv[1]);
	//int n2  = atoi(argv[2]);
	//MultiplyWithoutMult(n1, n2);

	//int numElements = atoi(argv[1]);
	//int k = atoi(argv[2]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//   elements[i] = atoi(argv[3+i]);
	//KthSmallestElement(elements, numElements,k);

	///Rain water trap
	//int numElements = atoi(argv[1]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//  elements[i] = atoi(argv[1+i+1]);
	//RainWaterTrap(elements,numElements);

	//Coin Sum
	//int numElements = atoi(argv[1]);
	//int targetSum = atoi(argv[2]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//   elements[i] = atoi(argv[3+i]);
	//CoinSum(elements,numElements, targetSum);

	//Linear time Median
	//int numElements = atoi(argv[1]);
	//int *elements = new int[numElements];
	//for (int i=0;i<numElements;i++)
	//  elements[i] = atoi(argv[1+i+1]);
	//MedianLinearTime(elements,numElements);


	unsigned int n1 = atoi(argv[1]);
	unsigned int n2  = atoi(argv[2]);
	MultiplyWithoutMultAdd(n1, n2);


	//delete [] elements;
}

unsigned int get_sum_no_add(unsigned int &a, unsigned int &b)
{
	unsigned int k = 1;
	unsigned int carry_in = 0;
	unsigned int sum = 0;
	unsigned int carry_out;

	while (k)
	{
		unsigned int ak = (a&k);
		unsigned int bk = (b&k);
		
		carry_out = (ak&bk) | (ak&carry_in) | (bk&carry_in);
		sum |= (ak^bk^carry_in);
		carry_in = carry_out << 1;

		k <<= 1;
	}
	return sum;
}

int MultiplyWithoutMultAdd(unsigned int &n1,unsigned int &n2)
{

	unsigned int k   = 1;
	unsigned int sum = 0;

	while (k)
	{
		if (n1&k)
			sum = get_sum_no_add(sum, n2);
		
		k <<= 1;
		n2 <<= 1;
	}

	cout << "Multiplication is :" << sum << endl;
	return 0;
}
void swap(int &A, int &B)
{
	int temp = A;
	A = B;
	B = temp;
}

int Partition(int *A, int left, int right, int pivot)
{
	int pivotValue = A[pivot];
	int partitionIndex = left;

	cout << "Inside partition function with (left, right, pivot):" << left << "," << right << "," << pivot << endl;

	swap(A[pivot], A[right]);
	for (int i = left; i < right; i++)
	{
		if (A[i] > pivotValue)
			swap(A[i], A[partitionIndex++]);
	}

	swap(A[right], A[partitionIndex]);
	return partitionIndex;
}

int MedianLinearTime(int *elements, int numElements)
{
	int left = 0;
	int right = numElements - 1;
	int k;

	if (numElements & 1)
		k = numElements / 2 ;
	else
		k = numElements / 2 -1;

	cout << "k is :" << k << endl;

	int median;
	while (left <= right)
	{
		default_random_engine gen;
		uniform_int_distribution<int> dis(left, right);
		int p = Partition(elements, left, right, dis(gen));
		if (p == k)
		{
			median = elements[p];
			break;
		}
		else if (p > k)
		{
			right = p - 1;
		}
		else
		{
			left = p + 1;
		}
	}

	cout << "Median value:" << median << endl;
	return 0;
}


int CoinSum(int *elem, int numElements, int targetSum)
{
	int i, j, x, y;
	int n = targetSum;
	int m = numElements;
	// We need n+1 rows as the table is constructed 
	
	// in bottom up manner using the base case 0
	// value case (n = 0)
	int table[100][100];

	// Fill the enteries for 0 value case (n = 0)
	for (i = 0; i<m; i++)
		table[0][i] = 1;

	// Fill rest of the table entries in bottom 
	// up manner  
	for (i = 1; i < n + 1; i++)
	{
		for (j = 0; j < m; j++)
		{
			// Count of solutions including S[j]
			x = (i - elem[j] >= 0) ? table[i - elem[j]][j] : 0;

			// Count of solutions excluding S[j]
			y = (j >= 1) ? table[i][j - 1] : 0;

			// total count
			table[i][j] = x + y;
		}
	}
	
	cout << "Number of possible ways:" << table[n][m - 1] << endl;
	for (i = 0; i < n + 1; i++)
	{
		for (j = 0; j < m; j++)
		{
			cout << table[i][j] << ",";
		}
		cout << endl;
	}

	return table[n][m - 1];

}

int RainWaterTrap(int *elem, int numElements)
{
	priority_queue<int> pq;
	for (int i = 0; i < numElements; i++)
		pq.push(elem[i]);
	
	int curr_max = pq.top();
	pq.pop();

	int curr_left = elem[0];
	int water_vol = 0;
	int start_loc = 0;
	int stop_loc  = 0;
	int max_level = 0;

	for (int i = 1; i < numElements; i++)
	{
		if (elem[i] < curr_left && i<numElements-1)
			stop_loc++;
		else
		{
			max_level = min(curr_left, elem[i]);
			for (int j = start_loc + 1; j <= stop_loc; j++)
				water_vol += max_level - elem[j];
			start_loc = i;
			curr_left = elem[i];
			stop_loc  = i;
		}
	}

	cout << "Water retained:" << water_vol << endl;
	return 0;
}
int KthSmallestElement(int *elem, int numElements, int k)
{
	priority_queue<int> pq;
	for (int i = 0; i < k; i++)
		pq.push(elem[i]);

	for (int i = k; i < numElements; i++)
	{
		if (elem[i] < pq.top())
		{
			pq.pop();
			pq.push(elem[i]);
		}
	}

	cout << "Kth largest number is:" << pq.top();
	return 0;
}

int MultiplyWithoutMult(int n1, int n2)
{
	int mult = 0;

	while (n1)
	{
		if ((n1 & 1) == 1)
			mult += n2;
		n1 >>= 1;//shift n1 to the right
		n2 <<= 1;//multiply n2 by2
	}
	cout << "Result of multiplication:" << mult << endl;
	return 0;
}

int TaskOptimization(int *elem, int numTasks, int maxTime)
{
	int temp;
	//first sort the element vector
    for(int i=0;i<numTasks-1;i++)
		for (int j = i + 1; j < numTasks; j++)
		{
			if (elem[j] < elem[i])
			{
				temp = elem[i];
				elem[i] = elem[j];
				elem[j] = temp;
			}
		}

	//now go through the tasks and accumulate
	int sum = 0;
	int count = 0;
	for (int i = 0; i < numTasks; i++)
	{
		if (sum+elem[i] <= maxTime)
		{ 
			sum += elem[i];
			count++;
		}
	}

	cout << "Number of tasks:" << count;
	return 0;
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



