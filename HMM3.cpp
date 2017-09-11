//Author Jonathan Rinnarv & Marcus Lignercrona

#include <iostream>
#include <vector>
#include <cmath>
class Matrix{
	int cols;
	int rows;
	std::vector<std::vector<long double>> m;
public:
	Matrix(){
		std::cin >> rows;
		std::cin >> cols;
		long double n;
		for(int r = 0; r < rows; ++r){
			std::vector<long double> row;
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
			std::vector<long double> row;
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

	long double& operator()(int row, int col){
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

long double calculateAlpha(Matrix& alpha, std::vector<int> & emission, Matrix& a, Matrix& b, Matrix& pi){
	//Alpha init
	for(int i = 0; i < alpha.height(); ++i){
		alpha(i,0) = pi(0,i)*b(i,emission[0]);
	}

	//Alpha loop
	for(int t = 1; t < emission.size(); ++t){
		for(int i = 0; i < alpha.height(); ++i){
			long double sum = 0;
			for(int j = 0; j < a.height(); ++j){
				sum += a(j,i)*alpha(j,t-1);
			}
			alpha(i,t) = b(i,emission[t])*sum;
		}
	}

	//Alpha close
	long double sum = 0;
	for(int i = 0; i < alpha.height(); ++i){
		sum+=alpha(i, emission.size()-1);
	}
	return sum;
}

void calculateBeta(Matrix& beta, std::vector<int> & emission, Matrix& a, Matrix& b){
	//Beta init
	for(int i = 0; i < beta.height(); ++i){
		beta(i,beta.width()-1) = 1;
	}
	
	//Beta loop
	for(int t = beta.width()-2; t >= 0; --t){
		for(int i = 0; i < beta.height(); ++i){
			long double sum = 0;			
			for(int j = 0; j < beta.height(); ++j){
				sum += beta(j,t+1)*b(j,emission[t+1])*a(i,j);
			}
			beta(i,t) = sum;
		}

	}
}

void calculateDiGamma(std::vector<Matrix>& diGamma, std::vector<int> & emission, Matrix& a, Matrix& b, Matrix& alpha, Matrix& beta, long double alphaSum){
	for(int t = 0; t < diGamma.size(); ++t){
		for(int i = 0; i < diGamma[0].height(); ++i){
			for(int j = 0; j < diGamma[0].width(); ++j){
				diGamma[t](i,j) = (alpha(i,t)*a(i,j)*b(j,emission[t+1])*beta(j,t+1))/alphaSum;	
			}
		}

	}
}

void getEmission(std::vector<int>& emission){
	int n;
	int e;
	std::cin >> n;
	for(int i = 0; i < n; ++i){
		std::cin >> e;
		emission.push_back(e);
	}
}

void calculateGamma(Matrix& gamma, std::vector<Matrix> & diGamma){
	for(int t = 0; t < diGamma.size(); ++t){
		for(int i = 0; i < diGamma[0].height(); ++i){
			long double sum = 0;
			for(int j = 0; j < diGamma[0].width(); ++j){
				sum+=diGamma[t](i,j);
			}
			gamma(i,t) = sum;
		}
	}
}

void updateValues(Matrix& a, Matrix& b, Matrix& pi, Matrix & gamma, std::vector<Matrix> & diGamma, std::vector<int> emission){
	//a
	for(int i = 0; i < a.height(); ++i){
		for(int j = 0; j < a.width(); ++j){
			long double diGammaSum = 0;
			long double gammaSum = 0;
			for(int t = 0; t < gamma.width(); ++t){
				diGammaSum+=diGamma[t](i,j);
				gammaSum += gamma(i,t);
			}
			a(i,j)=diGammaSum/gammaSum;
		}
	}

	//b
	for(int j = 0; j < b.height(); ++j){
		for(int k = 0; k < b.width(); ++k){
			long double numerator = 0;
			long double denom = 0;
			for(int t = 0; t < gamma.width(); ++t){
				if(emission[t] == k)
					numerator += gamma(j,t);
				denom += gamma(j,t);
			}
			b(j,k)=numerator/denom;
		}
	}
	
	//pi
	for(int i = 0; i < pi.width(); ++i)
		pi(0,i) = gamma(i,0);
}

bool continueEpsilon(Matrix& alpha, long double & oldLogProb){
	long double logProb = 0;
	for(int t = 0; t < alpha.width(); ++t){
		long double sum = 0;
		for(int i = 0; i < alpha.height(); ++i){
			sum+=alpha(i,t);
		}
		logProb += std::log((long double)1/sum);
	}
	//std::cout << logProb-oldLogProb << std::endl;
	if(std::abs(logProb-oldLogProb) < 0.01)
		return true;
	oldLogProb = logProb;
	return false;
}

main(){
	Matrix a = Matrix();
	Matrix b = Matrix();
	Matrix pi = Matrix();

	std::vector<int> emission;
	getEmission(emission);
	long double oldLogProb = 0;
	
	for(int i = 0; i < 30; ++i){
		Matrix alpha = Matrix(a.height(), emission.size());
		long double alphaSum = calculateAlpha(alpha, emission, a, b, pi);

		Matrix beta = Matrix(a.height(), emission.size());
		calculateBeta(beta, emission, a, b);
	
		std::vector<Matrix> diGamma;
		Matrix gammaTemp = Matrix(a.height(), a.width());
		for(int t = 0; t < emission.size()-1; ++t)		 
			diGamma.push_back(gammaTemp);

	
		calculateDiGamma(diGamma, emission, a, b, alpha, beta, alphaSum);

		Matrix gamma = Matrix(a.height(), emission.size()-1);
		calculateGamma(gamma,diGamma);

		updateValues(a,b,pi,gamma,diGamma, emission);
		if(continueEpsilon(alpha, oldLogProb)){
			//std::cout << i << " " << std::endl;
			break;
		}
		
	}
	
	std::cout << a.height() << " " << a.width() << " ";
	for(int i = 0; i < a.height(); ++i){
		for(int j = 0; j < a.width(); ++j){
			std::cout << a(i,j) << " ";
		}
	}

	std::cout << std::endl;
	std::cout << b.height() << " " << b.width() << " ";
	for(int i = 0; i < b.height(); ++i){
		for(int j = 0; j < b.width(); ++j){
			std::cout << b(i,j) << " ";
		}
	}
	std::cout << std::endl;
		
}
