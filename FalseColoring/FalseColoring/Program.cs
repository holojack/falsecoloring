using System;
using System.Diagnostics; // For the stopwatch
using System.Threading.Tasks; // For the Parallel

namespace FalseColoring
{
    class Program
    {
        static void Main(string[] args)
        {
            int width = 1000;
            int length = 5000;

            // Declare the 2D array
            int[,] grid = new int[width, length];

            // Fill the 2D array with the starting values
            grid = SetupGrid(grid);

            // Program is done, make sure to the let user see the results
            System.Console.WriteLine("Press Any Key to Exit");
            Console.ReadKey();

        }

        static int[,] SetupGrid(int[,] grid)
        {
            // Loop over the columns
            for (int col = 0; col < grid.GetLength(0); col++)
            {
                // Loop over the rows
                for (int row = 0; row < grid.GetLength(1); row++)
                {
                    if (col == 0)
                    {
                        // The far left column of the grid is 2000 degrees celsius by default
                        grid[col, row] = 2000;
                    }
                    else
                    {
                        // Otherwise, start with zero degrees celsius
                        grid[col, row] = 0;
                    }
                }
            }

            return grid;
        }
    }
}
