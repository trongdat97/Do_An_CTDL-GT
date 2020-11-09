#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
using namespace std;

const int N = 10000;
int arrRandom[N], arrSorted[N], arrReversed[N];

// A utility function to swap two elements 
void swap(int &x, int &y)  
{  
    int temp = x;  
    x = y;  
    y = temp;  
}  

//===============Selection Sort===============// 
float selectionSort(int arr[], int n)  
{  
	clock_t begin = clock();   
    int i, j, min_idx;  
  
    // One by one move boundary of unsorted subarray  
    for (i = 0; i < n-1; i++)  
    {  
        // Find the minimum element in unsorted array  
        min_idx = i;  
        for (j = i+1; j < n; j++)  
        if (arr[j] < arr[min_idx])  
            min_idx = j;  
  
        // Swap the found minimum element with the first element  
        swap(arr[min_idx], arr[i]);  
    }  
    
    clock_t end = clock(); 
    return (float)(end-begin)/CLOCKS_PER_SEC;
}  
//============================================//

//===============Insertion Sort===============//
float insertionSort(int arr[], int n)
{
	clock_t begin = clock(); 
   	int i, key, j;  
    for (i = 1; i < n; i++) 
    {  
        key = arr[i];  
        j = i - 1;  
  
        /* Move elements of arr[0..i-1], that are  
        greater than key, to one position ahead  
        of their current position */
        while (j >= 0 && arr[j] > key) 
        {  
            arr[j + 1] = arr[j];  
            j = j - 1;  
        }  
        arr[j + 1] = key;  
    }  
   	clock_t end = clock(); 
    return (float)(end-begin)/CLOCKS_PER_SEC;
}
//============================================//

//===============Bubble Sort===============//
float bubbleSort(int arr[], int n)  
{  
	clock_t begin = clock(); 
    int i, j;  
    for (i = 0; i < n-1; i++)      
      
    // Last i elements are already in place  
    for (j = 0; j < n-i-1; j++)  
        if (arr[j] > arr[j+1])  
            swap(arr[j], arr[j+1]);  
    clock_t end = clock(); 
    return (float)(end-begin)/CLOCKS_PER_SEC;
}  
//============================================//

//===============Quick Sort===================//

/* The main function that implements QuickSort  
arr[] --> Array to be sorted,  
low --> Starting index,  
high --> Ending index */
float quickSort(int arr[], int l, int r)  
{  
	clock_t begin = clock(); 
    // If the first index less or equal than the last index
	if (l <= r)
	{
		// Create a Key/Pivot Element
		int key = arr[(l+r)/2];
 
		// Create temp Variables to loop through array
		int i = l;
		int j = r;
 
		while (i <= j)
		{
			while (arr[i] < key)
				i++;
			while (arr[j] > key)
				j--;
 
			if (i <= j)
			{
				swap(arr[i], arr[j]);
				i++;
				j--;
			}
		}
 
		// Recursion to the smaller partition in the array after sorted above
		if (l < j)
			quickSort(arr, l, j);
		if (r > i)
			quickSort(arr, i, r);
	}
    clock_t end = clock(); 
    return (float)(end-begin)/CLOCKS_PER_SEC;
}  
//============================================//

//=================Merge Sort=================//
// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge(int arr[], int l, int m, int r) 
{ 
    int i, j, k; 
    int n1 = m - l + 1; 
    int n2 =  r - m; 
  
    /* create temp arrays */
    int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++) 
        L[i] = arr[l + i]; 
    for (j = 0; j < n2; j++) 
        R[j] = arr[m + 1+ j]; 
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 
    while (i < n1 && j < n2) 
        arr[k++] = (L[i] <= R[j] ? L[i++] : R[j++]);
        
    /* Copy the remaining elements of L[], if there are any */
    while (i < n1) 
    { 
        arr[k++] = L[i++]; 
         
    } 
  
    /* Copy the remaining elements of R[], if there are any */
    while (j < n2) 
    { 
        arr[k++] = R[j++]; 
    } 
} 

/* l is for left index and r is right index of the sub-array of arr to be sorted */
float mergeSort(int arr[], int l, int r) 
{ 
    clock_t begin = clock(); 
	if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        int m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort(arr, l, m); 
        mergeSort(arr, m+1, r); 
  
        merge(arr, l, m, r); 
    } 
    clock_t end = clock(); 
    return (float)(end-begin)/CLOCKS_PER_SEC;
} 
//============================================//

//================Shell Sort==================//
float shellSort(int arr[], int n) 
{ 
	clock_t begin = clock(); 
    for (int interval = n/2; interval > 0; interval /= 2) 
    { 
        // Do a intervalped insertion sort for this interval size. 
        // The first interval elements a[0..interval-1] are already in intervalped order 
        // keep adding one more element until the entire array is 
        // interval sorted  
        for (int i = interval; i < n; i += 1) 
        { 
            // add a[i] to the elements that have been interval sorted 
            // save a[i] in temp and make a hole at position i 
            int temp = arr[i]; 
  
            // shift earlier interval-sorted elements up until the correct  
            // location for a[i] is found 
            int j;             
            for (j = i; j >= interval && arr[j - interval] > temp; j -= interval) 
                arr[j] = arr[j - interval]; 
              
            //  put temp (the original a[i]) in its correct location 
            arr[j] = temp; 
        } 
    } 
    clock_t end = clock(); 
    return (float)(end-begin)/CLOCKS_PER_SEC;
} 
//============================================//

//=================Radix Sort=================//
// A utility function to get maximum value in arr[] 
int getMax(int arr[], int n) 
{ 
    int mx = arr[0]; 
    for (int i = 1; i < n; i++) 
        if (arr[i] > mx) 
            mx = arr[i]; 
    return mx; 
} 

// A function to do counting sort of arr[] according to 
// The digit represented by exp.
void countSort(int arr[], int n, int exp) 
{ 
    int output[n]; // output array 
    int i, count[10] = {0}; 
  
    // Store count of occurrences in count[] 
    for (i = 0; i < n; i++) 
        count[ (arr[i]/exp)%10 ]++; 
  
    // Change count[i] so that count[i] now contains actual 
    //  position of this digit in output[] 
    for (i = 1; i < 10; i++) 
        count[i] += count[i - 1]; 
  
    // Build the output array 
    for (i = n - 1; i >= 0; i--) 
    { 
        output[count[ (arr[i]/exp)%10 ] - 1] = arr[i]; 
        count[ (arr[i]/exp)%10 ]--; 
    } 
  
    // Copy the output array to arr[], so that arr[] now 
    // contains sorted numbers according to current digit 
    for (i = 0; i < n; i++) 
        arr[i] = output[i]; 
} 

// The main function to that sorts arr[] of size n using  
// Radix Sort 
float radixSort(int arr[], int n) 
{ 
	clock_t begin = clock(); 
    // Find the Nimum number to know number of digits 
    int m = getMax(arr, n); 
  
    // Do counting sort for every digit. Note that instead 
    // of passing digit number, exp is passed. exp is 10^i 
    // where i is current digit number 
    for (int exp = 1; m/exp > 0; exp *= 10) 
        countSort(arr, n, exp); 
    clock_t end = clock(); 
    return (float)(end-begin)/CLOCKS_PER_SEC;
} 
//============================================//

//=================Heap sort=================//
// To heapify a subtree rooted with node i which is 
// an index in arr[]. n is size of heap 
void heapify(int arr[], int n, int i) 
{ 
    int largest = i; // Initialize largest as root 
    int l = 2*i + 1; // left = 2*i + 1 
    int r = 2*i + 2; // right = 2*i + 2 
  
    // If left child is larger than root 
    if (l < n && arr[l] > arr[largest]) 
        largest = l; 
  
    // If right child is larger than largest so far 
    if (r < n && arr[r] > arr[largest]) 
        largest = r; 
  
    // If largest is not root 
    if (largest != i) 
    { 
        swap(arr[i], arr[largest]); 
  
        // Recursively heapify the affected sub-tree 
        heapify(arr, n, largest); 
    } 
} 

//Main function to do heap sort 
float heapSort(int arr[], int n) 
{ 
	clock_t begin = clock(); 
    // Build heap (rearrange array) 
    for (int i = n / 2 - 1; i >= 0; i--) 
        heapify(arr, n, i); 
  
    // One by one extract an element from heap 
    for (int i=n-1; i>=0; i--) 
    { 
        // Move current root to end 
        swap(arr[0], arr[i]); 
  
        // call N heapify on the reduced heap 
        heapify(arr, i, 0); 
    } 
    clock_t end = clock(); 
    return (float)(end-begin)/CLOCKS_PER_SEC;
} 
//============================================//

// A utility function to copy array
void copyArray(int arr1[], int arr2[], int n)
{
	for (int i=0; i<n; i++) arr1[i] = arr2[i];
}

//This function to read input file
void readFile(){
	//Input the file was created by the CreateInput.cpp
	ifstream infile;
	infile.open("input.txt");
	for (int i = 0; i<N; i++ ) infile>>arrRandom[i];
	
	int temp[N];
	//Copy Array arrRandom to Array temp
	copyArray(temp, arrRandom, N);
	//Sort Array temp
	quickSort(temp, 0, N-1);
	//Copy Array tempp to Array arrSorted
	copyArray(arrSorted, temp, N);
	//Reverse arrSorted and assign it to arrReversed
	for (int i=0; i<N; i++) arrReversed[N-i-1] = arrSorted[i];
    infile.close();
}

int main(){
	
	int arr[N];
	//Call the read file function
	readFile();

	ofstream outfile;
	outfile.open("output.txt");
	
	//-----------------------------Array with Random values----------------------------------
	outfile << setfill('=');
	outfile << setw(18) << "=";
	outfile<<"Random";
	outfile << setw(18) << "="<<endl<<endl;
	
	//Write the result of SelectionSort to the output file
	outfile << setfill(' ');
	copyArray(arr, arrRandom, N);
	outfile << setw(20) << left << " - Selection Sort" << ": ";
	outfile << setw(8) << left << setprecision(6) << fixed << selectionSort(arr, N);
	outfile << setw(10) << right << " seconds"<<endl;

	//Write the result of InsertionSort to the output file
	copyArray(arr, arrRandom, N);
	outfile << setw(20) << left << " - Insertion Sort" << ": ";
	outfile << setw(8) << left << insertionSort(arr, N);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of BubbleSort to the output file
	copyArray(arr, arrRandom, N);
	outfile << setw(20) << left << " - Bubble Sort" << ": ";
	outfile << setw(8) << left << bubbleSort(arr, N);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of QuickSort to the output file
	copyArray(arr, arrRandom, N);
	outfile << setw(20) << left << " - Quick Sort" << ": ";
	outfile << setw(8) << left << quickSort(arr, 0, N-1);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of MergeSort to the output file	
	copyArray(arr, arrRandom, N);
	outfile << setw(20) << left << " - Merge Sort" << ": ";
	outfile << setw(8) << left << mergeSort(arr, 0, N-1);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of ShellSort to the output file	
	copyArray(arr, arrRandom, N);
	outfile << setw(20) << left << " - Shell Sort" << ": ";
	outfile << setw(8) << left << shellSort(arr, N);
	outfile << setw(10) << right << " seconds" << " " << endl;

	//Write the result of RadixSort to the output file
	copyArray(arr, arrRandom, N);
	outfile << setw(20) << left << " - Radix Sort" << ": ";
	outfile << setw(8) << left << radixSort(arr, N);
	outfile << setw(10) << right << " seconds" << " " << endl;
	
	//Write the result of HeapSort to the output file	
	copyArray(arr, arrRandom, N);
	outfile << setw(20) << left << " - Heap Sort" << ": ";
	outfile << setw(8) << left << heapSort(arr, N);
	outfile << setw(10) << right << " seconds" << " " << endl<<endl;
	
	outfile << setfill('=');		
	outfile << setw(42) << "=" << endl<<endl;
	
	//-----------------------------------Array with Sorted values------------------------------
	
	outfile << setw(18) << "=";
	outfile<<"Sorted";
	outfile << setw(18) << "="<<endl<<endl;
	
	//Write the result of SelectionSort to the output file
	outfile << setfill(' ');
	copyArray(arr, arrSorted, N);
	outfile << setw(20) << left << " - Selection Sort" << ": ";
	outfile << setw(8) << left << selectionSort(arr, N);
	outfile << setw(10) << right << " seconds"<<endl;

	//Write the result of Insertion to the output file
	copyArray(arr, arrSorted, N);
	outfile << setw(20) << left << " - Insertion Sort" << ": ";
	outfile << setw(8) << left << insertionSort(arr, N);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of BubbleSort to the output file
	copyArray(arr, arrSorted, N);
	outfile << setw(20) << left << " - Bubble Sort" << ": ";
	outfile << setw(8) << left << bubbleSort(arr, N);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of QuickSort to the output file
	copyArray(arr, arrSorted, N);
	outfile << setw(20) << left << " - Quick Sort" << ": ";
	outfile << setw(8) << left << quickSort(arr, 0, N-1);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of MergeSort to the output file	
	copyArray(arr, arrSorted, N);
	outfile << setw(20) << left << " - Merge Sort" << ": ";
	outfile << setw(8) << left << mergeSort(arr, 0, N-1);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of ShellSort to the output file
	copyArray(arr, arrSorted, N);
	outfile << setw(20) << left << " - Shell Sort" << ": ";
	outfile << setw(8) << left << shellSort(arr, N);
	outfile << setw(10) << right << " seconds" << " " << endl;

	//Write the result of RadixSort to the output file
	copyArray(arr, arrSorted, N);
	outfile << setw(20) << left << " - Radix Sort" << ": ";
	outfile << setw(8) << left << radixSort(arr, N);
	outfile << setw(10) << right << " seconds" << " " << endl;
	
	//Write the result of HeapSort to the output file
	copyArray(arr, arrSorted, N);
	outfile << setw(20) << left << " - Heap Sort" << ": ";
	outfile << setw(8) << left << heapSort(arr, N);
	outfile << setw(10) << right << " seconds" << " " << endl<<endl;
	
	outfile << setfill('=');		
	outfile << setw(42) << "=" << endl<<endl;
	
	//----------------------------Array with Reversed values--------------------------------

	outfile << setw(17) << "=";
	outfile<<"Reversed";
	outfile << setw(17) << "="<<endl<<endl;
	
	//Write the result of SelectionSort to the output file
	outfile << setfill(' ');
	copyArray(arr, arrReversed, N);
	outfile << setw(20) << left << " - Selection Sort" << ": ";
	outfile << setw(8) << left << selectionSort(arr, N);
	outfile << setw(10) << right << " seconds"<<endl;

	//Write the result of Insertion to the output file	
	copyArray(arr, arrReversed, N);
	outfile << setw(20) << left << " - Insertion Sort" << ": ";
	outfile << setw(8) << left << insertionSort(arr, N);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of BubbleSort to the output file
	copyArray(arr, arrReversed, N);
	outfile << setw(20) << left << " - Bubble Sort" << ": ";
	outfile << setw(8) << left << bubbleSort(arr, N);
	outfile << setw(10) << right << " seconds" << endl;
	
    //Write the result of QuickSort to the output file
	copyArray(arr, arrReversed, N);
	outfile << setw(20) << left << " - Quick Sort" << ": ";
	outfile << setw(8) << left << quickSort(arr, 0, N-1);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of MergeSort to the output file	
	copyArray(arr, arrReversed, N);
	outfile << setw(20) << left << " - Merge Sort" << ": ";
	outfile << setw(8) << left << mergeSort(arr, 0, N-1);
	outfile << setw(10) << right << " seconds" << endl;
	
	//Write the result of ShellSort to the output file
	copyArray(arr, arrReversed, N);
	outfile << setw(20) << left << " - Shell Sort" << ": ";
	outfile << setw(8) << left << shellSort(arr, N);
	outfile << setw(10) << right << " seconds" << " " << endl;

	//Write the result of RadixSort to the output file
	copyArray(arr, arrReversed, N);
	outfile << setw(20) << left << " - Radix Sort" << ": ";
	outfile << setw(8) << left << radixSort(arr, N);
	outfile << setw(10) << right << " seconds" << " " << endl;

	//Write the result of HeapSort to the output file
	copyArray(arr, arrReversed, N);
	outfile << setw(20) << left << " - Heap Sort" << ": ";
	outfile << setw(8) << left << heapSort(arr, N);
	outfile << setw(10) << right << " seconds" << " " << endl<<endl;
	
	outfile << setfill('=');		
	outfile << setw(42) << "=" << endl;
	
	cout<<"Toan bo ket qua da duoc luu vao file output.txt!";
	outfile.close();
	
    return 0;
}
