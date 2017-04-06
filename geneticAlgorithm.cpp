#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <vector>
#include <float.h>
#include <iostream>
#include <bitset>
#include <exception>
#include <math.h>

using function = double (*)(double);

int randrange(int min, int max){
	return min + (rand() % (int)(max - min + 1));
}

double my_func(double x){
	return x*sin(10*M_PI * x) + 1;
}

class Population{
	private:
		int chromosome_size = 22;
		double p_mut;
		double p_cross;
		double step;
		int beg, end;
		int pop_size;
		std::vector<int> population;
		std::vector<double> weights;
		int best_value;
		double best_result;
		function f;
		void CreatePopulation();
		void FindIfBest();
		void Roulette();
		void Crossover();
		void Mutation();
		int BitFlipMutation(int value);
		int SinglePointCrossover(int parent_first, int parent_second);
		double StepToValue(int i);
		int ValueToStep(double x);
	public:
		Population(const int n, const int a, const int b, function _f);
		virtual ~Population();
		std::vector<int> GetPopulation() const;
		void Init(const int number_of_iterations);
        double getBestResult();
        double getBestValue();
        void setP_Cross(double _p_cross);
        void setP_Mut(double _p_mut);
};

Population::Population(const int n, const int a, const int b, function _f){
			p_mut = 0.2;
			p_cross = 0.8;
			f = _f;
			best_result = -DBL_MAX;
			pop_size = n;
			beg = a;
			end = b;
			CreatePopulation();
}

std::vector<int> Population::GetPopulation() const{
	return population;
}


void Population::setP_Cross(const double _p_cross){
	p_cross = _p_cross;
}

void Population::setP_Mut(const double _p_mut){
	p_mut = _p_mut;
}

void Population::CreatePopulation(){
	int steps = pow(2,chromosome_size);
	step = (double)(end-beg)/(double)steps;
	int step_value;
	for(int i = 0; i < pop_size; i++){
		step_value = randrange(0,steps);
		population.push_back(step_value);
		weights.push_back(f(StepToValue(step_value)));
	}
	FindIfBest();
}

void Population::FindIfBest(){
	for(std::vector<int>::size_type i = 0; i != population.size() && i != weights.size(); i++){
		//std::cout<<best_result<<" : "<<weights[i]<<"\n";
		if(best_result < weights[i])
		{
			best_value = population[i];
			best_result = weights[i];
		}
	}
}

void Population::Roulette(){
	double sum_of_weights = 0.0;
	for(auto const& value: weights)
		sum_of_weights += value;
	std::vector<double> normalized_weights;

	for(auto const& value: weights)
		normalized_weights.push_back(value/sum_of_weights);

	for(std::vector<int>::size_type i = 0; i != normalized_weights.size(); i++){
			if(i > 0)
				normalized_weights[i] += normalized_weights[i-1];
	}
	double r;
	int u = 0;
	std::vector<int> temp_population;
	for(std::vector<int>::size_type i = 0; i != population.size(); i++){
			r = ((double) rand() / (RAND_MAX));
			while(r >= normalized_weights[u])
				u++;
			temp_population.push_back(population[u]);
			u = 0;
	}
	population = temp_population;
}

void Population::Mutation(){
	std::vector<int> to_mutate;
	double mutate_size = p_mut*chromosome_size;
	int random_index;
	for(unsigned int i = 0; i < mutate_size; i++){
		random_index = randrange(0,population.size()-1);
		to_mutate.push_back(population[random_index]);
		population.erase(population.begin()+random_index);
	}
	for(std::vector<int>::size_type i = 0; i != to_mutate.size(); i++){
        for(unsigned int u = 0; u < mutate_size; u++)
            to_mutate[i] = BitFlipMutation(to_mutate[i]);
        population.push_back(to_mutate[i]);
	}
}

int Population::BitFlipMutation(int value){
	int random_index = randrange(0,chromosome_size-1);
	int mask = 1<<random_index;
	int result = value ^ mask;
	return result;
}


void Population::Crossover(){
	std::vector<int > to_crossover;
	double crossover_size = p_cross*pop_size;
	int random_index;
	for(int i = 0; i < crossover_size; i++){
		random_index = randrange(0,population.size()-1);
		to_crossover.push_back(population[random_index]);
		population.erase(population.begin()+random_index);
	}

	if(to_crossover.size()%2 != 0){
		random_index = randrange(0,population.size()-1);
		to_crossover.push_back(population[random_index]);
		population.erase(population.begin()+random_index);
	}

	int parent_first, parent_second;
	int child_first, child_second;
	std::vector<int> after_crossover;
    unsigned int limit = to_crossover.size()/2;
	for(unsigned int i = 0; i < limit; i++){
		//first parent
		random_index = randrange(0,to_crossover.size()-1);
		parent_first = to_crossover[random_index];
		to_crossover.erase(to_crossover.begin()+random_index);
		//second parent
		random_index = randrange(0,to_crossover.size()-1);
		parent_second = to_crossover[random_index];
		to_crossover.erase(to_crossover.begin()+random_index);


		child_first = SinglePointCrossover(parent_first, parent_second);
		child_second = SinglePointCrossover(parent_second, parent_first);

		after_crossover.push_back(child_first);
		after_crossover.push_back(child_second);
	}

	population.insert(population.begin(), after_crossover.begin(), after_crossover.end());
}

double Population::StepToValue(int i){
	return beg + step*i;
}

int Population::ValueToStep(double x){
	int to_return = (x - beg)/step;
	return to_return;
}

int Population::SinglePointCrossover(int parent_first, int parent_second){
	int x = parent_first;
	int y = parent_second;
	int random_index = randrange(0,chromosome_size);

	int tmp_x = x >> random_index;
	tmp_x = tmp_x << random_index;
	int mask = (1 << random_index) - 1;
	int tmp_y = y & mask;
	int result = tmp_x | tmp_y;

	return result;
}

void Population::Init(const int number_of_iterations){

	for(int i = 0; i < number_of_iterations; i++){
		try{
            //selection
			Roulette();
			//crossing
			Crossover();
			//mutate
			Mutation();
			//judge
			for(std::vector<int>::size_type i = 0; i != population.size() && i != weights.size(); i++)
                weights[i] = f(StepToValue(population[i]));
			FindIfBest();
		}catch(std::exception& e){
			std::cout<<e.what()<<'\n';
		}
	}
}

Population::~Population(){

}

double Population::getBestValue() {
    return StepToValue(best_value);
}


double Population::getBestResult() {
    return best_result;
}

double test(function f, int beg, int end){
    double max = -DBL_MAX;
    unsigned int steps = pow(2,22);
    double step = (double )(end - beg)/(double )steps;
    for(unsigned int i = 0; i < steps; i++){
        if(max < f(beg + step*i))
            max = f(beg + step*i);
    }
    return  max;
}


int main(int argc, char const *argv[])
{
	srand( time( NULL ) );
	//temp.push_back(randrange(0,pow(2,22)));
    std::clock_t  start;
    double duration_genetic, duration_analitic;
	Population geteticTest(20, -1, 2, my_func);
    double test_result;

    start = std::clock();
	geteticTest.Init(100);
    duration_genetic = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    start = std::clock();
    test_result = test(my_func,-1,2);
    duration_analitic = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"\nGenetic Result: "<<geteticTest.getBestResult()<<"\nGenetic Time: "<<duration_genetic<<'\n';
    std::cout<<"\nAnalitic Result: "<<test_result<<"\nAnalitic Time: "<<duration_analitic<<'\n';

	/* code */
	return 0;
}