//Author Jonathan Rinnarv & Marcus Lignercrona

#include <iostream>
#include <vector>
#include <algorithm>

class Matrix{
	int cols;
	int rows;
	std::vector<std::vector<double>> m;
public:
	Matrix(){
		std::cin >> rows;
		std::cin >> cols;
		double n;
		for(int r = 0; r < rows; ++r){
			std::vector<double> row;
			for(int c = 0; c < cols; ++c){
				std::cin >> n;
				row.push_back(n);
			}
			m.push_back(row);
		}
	}

	Matrix(int row, int col){
		rows = row;
		cols = col;
		for(int r = 0; r < rows; ++r){
			std::vector<double> row;
			for(int c = 0; c < cols; ++c){
				row.push_back(-1);
			}
			m.push_back(row);
		}
	}

	int width(){
		return cols;
	}

	int height(){
		return rows;
	}

	double& operator()(int row, int col){
		return m[row][col];
	}

	void print(){
		std::cout << "--------------" << std::endl;
		for(int r = 0; r < rows; ++r){
			for(int c = 0; c < cols; ++c){
				std::cout << m[r][c] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "--------------" << std::endl;
	}

};


main(){
	Matrix a = Matrix();
	Matrix b = Matrix();
	Matrix pi = Matrix();

	std::vector<int> emission;
	int n;
	int e;
	std::cin >> n;
	for(int i = 0; i < n; ++i){
		std::cin >> e;
		emission.push_back(e);
	}

	Matrix delta = Matrix(a.height(), emission.size());
	Matrix deltaIndex = Matrix(a.height(), emission.size());


	//Alpha init
	for(int i = 0; i < delta.height(); ++i){
		delta(i,0) = pi(0,i)*b(i,emission[0]);
	}

	//delta loop
	for(int t = 1; t < emission.size(); ++t){
		for(int i = 0; i < delta.height(); ++i){
			double currentMax = -1;
			int currentArgMax = -1;
			double mightBeNewMax = -1;
			for(int j = 0; j < a.height(); ++j){
				mightBeNewMax = std::max(currentMax,a(j,i)*delta(j,t-1)*b(i,emission[t]));
				if(mightBeNewMax > currentMax){
					currentMax = mightBeNewMax;
					currentArgMax = j;
				}
			}
			delta(i,t) = currentMax;
			deltaIndex(i,t) = currentArgMax;
		}
	}

	//Alpha close
	std::vector<int> states;

	double largestProbability = -1;
	double mightBelargestProbability = -1;
	int index = -1;

	//set index of most probable last state
	for(int j = 0; j < a.height(); ++j){
		mightBelargestProbability = std::max(largestProbability,delta(j,emission.size()-1));
		if(mightBelargestProbability > largestProbability){
			largestProbability = mightBelargestProbability;
			index = j;
		}
	}
	//store most probable last state
	states.push_back(index);
	int prev = index;
	//go down the path
	for(int t = emission.size()-1;t>0;--t){
		prev = deltaIndex(prev,t);
		states.push_back(prev);
	}

	//output most probable path
	while(!states.empty()){
		std::cout << states.back() << " ";
		states.pop_back();
	}
}
