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

	std::vector<double> newPi;
	for(int j = 0; j < a.width(); ++j){
		double res = 0;	
		for(int i = 0; i < pi.width(); ++i){			
			res+=a(i,j)*pi(0,i);
		}
		newPi.push_back(res);
	}
	
	std::cout << "1 " << b.width();
	for(int k = 0; k < b.width(); ++k){
		double res = 0;
		for(int i = 0; i < newPi.size(); ++i){
			res += newPi[i]*b(i,k);
		}
		std::cout << " " << res;
	}
	std::cout << std::endl;
}
