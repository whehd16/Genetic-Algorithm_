#include "ga.h"
#include "simulation.h"

int main()
{
    srand(time(0));
    set_input_parameters();
    Population pop;
}

/*
int main()
{
    for (int i = 0; i < ITERATION; i++)
    {
        srand(time(0));
        set_input_parameters();
        Population pop;
        pop.run();
        pop.statistics_info_calc();
        cout << "<------------" << i << "th iteration------------>" << endl;
        /*cout << "assignment: " << endl;

		for (int j = 0; j < K; j++)
		{
			cout << assignment[j] << endl;
		}
		cout << endl;
        // pop.print_population(MODE_DETAIL);
        pop.print_best_assignment();
    }
}
*/
