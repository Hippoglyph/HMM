//Author Jonathan Rinnarv & Marcus Lignercrona

#include <iostream>
#include <vector>

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

	Matrix alpha = Matrix(a.height(), emission.size());

	//Alpha init
	for(int i = 0; i < alpha.height(); ++i){
		alpha(i,0) = pi(0,i)*b(i,emission[0]);
	}

	//Alpha loop
	for(int t = 1; t < emission.size(); ++t){
		for(int i = 0; i < alpha.height(); ++i){
			double sum = 0;
			for(int j = 0; j < a.height(); ++j){
				sum += a(j,i)*alpha(j,t-1);
			}
			alpha(i,t) = b(i,emission[t])*sum;
		}
	}

	//Alpha close
	double sum = 0;
	for(int i = 0; i < alpha.height(); ++i){
		sum+=alpha(i, emission.size()-1);
	}
	std::cout << sum << std::endl;
}
