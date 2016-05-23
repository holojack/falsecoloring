using System;
using System.Diagnostics; // For the stopwatch
using System.Threading.Tasks; // For the Parallel
using System.Linq; // For the 2D array comparison

namespace FalseColoring
{
    class Program
    {

        static void Main(string[] args)
        {
            // Set up the grid dimensions
            // Parallel will be slower for smaller grids
            int width = 12;
            int length = 12;

            // Declare the 2D array
            int[,] grid = new int[length, width];

            // Keep track of the run time
            Stopwatch sw = new Stopwatch();

            // Fill the 2D array with the starting values
            SetupGrid(grid);

            // Solve sequentially
            sw.Start();
            SeqGaussSeidel(grid);
            sw.Stop();
            Console.WriteLine("Elapsed time taken sequentially: {0}", sw.ElapsedMilliseconds);

            // Print the results
            for (int row = 0; row < width; row++)
            {
                for (int col = 0; col < length; col++)
                {
                    if (col == length - 1)
                    {
                        System.Console.WriteLine(grid[col, row]);
                    }
                    else
                    {
                        System.Console.Write(grid[col, row] + "-");
                    }
                }
            }

            // Fill the 2D array with the starting values
            SetupGrid(grid);

            // Solve in parallel
            sw.Reset();
            sw.Start();
            ParGaussSeidel(grid);
            sw.Stop();
            Console.WriteLine("Elapsed time taken parallel: {0}", sw.ElapsedMilliseconds);

            // Print the results
            for (int row = 0; row < width; row++)
            {
                for (int col = 0; col < length; col++)
                {
                    if (col == length - 1)
                    {
                        System.Console.WriteLine(grid[col, row]);
                    }
                    else
                    {
                        System.Console.Write(grid[col, row] + "-");
                    }
                }
            }

            // Program is done, make sure to the let user see the results
            System.Console.WriteLine("Press Any Key to Exit");
            Console.ReadKey();

        }

        /*
         * SeqTempSolver
         * 
         * Sets the starting temperature for each of the nodes assuming that the leftmost column
         * has a temperature of 2000 degrees celsius and everything else has a temperature of
         * 0 degrees celsius
         * 
         * Parameters:
         * int[,] grid - The input grid of all the nodes
         */
        static void SetupGrid(int[,] grid)
        {
            // Loop over the columns
            for (int col = 0; col < grid.GetLength(0); col++)
            {
                // Loop over the rows
                for (int row = 0; row < grid.GetLength(1); row++)
                {
                    if (col == 0)
                    {
                        // Have to initialize the left side to 2000 degrees
                        grid[0, row] = 2000;
                    }
                    else
                    {
                        // Otherwise, assume the temperature is 0 degrees
                        grid[col, row] = 0;
                    }

                    // This leaves out the corners, but those technically don't matter for the calculations
                }
            }
        }

        /*
         * SeqTempSolver
         * 
         * Uses a sequential version of the Gauss-Seidell algorithm to find the temperature of each
         * section of a 2D metal plate that is submerged in a solution with a temperature of 0 degrees
         * celsius on each of the sides with the exception of the left which has a temperature of
         * 2000 degrees celsius
         * 
         * Parameters:
         * int[,] grid - The input grid of all the nodes that temperatures will be determined for
         */
        static void SeqGaussSeidel(int[,] grid)
        {
            // Keep calculating the temperatures until they stabilize

            // Previous Grid
            // Used to make sure the final temperature has been found
            int[,] prevGrid = new int[grid.GetLength(0), grid.GetLength(1)];

            // Use LINQ to compare the grids (http://stackoverflow.com/questions/12446770/how-to-compare-multidimensional-arrays-in-c-sharp)
            // Only stop when the previous is the same as the current
            while (!(prevGrid.Rank == grid.Rank &&
                     Enumerable.Range(0, prevGrid.Rank).All(dimension => prevGrid.GetLength(dimension) == grid.GetLength(dimension)) &&
                     prevGrid.Cast<int>().SequenceEqual(grid.Cast<int>())))
            {
                // Cannot just assign it straight across because that passes the memory space not the values
                prevGrid = (int[,]) grid.Clone();

                // Calculate the temperatures
                SeqTempSolver(grid);
            }
        }

        /*
         * SeqTempSolver
         * 
         * Finds the temperature for each of the nodes
         * 
         * Parameters:
         * int[,] grid - The input grid of all the nodes that temperatures will be calculated for
         */
        static void SeqTempSolver(int[,] grid)
        {
            // Loop over the rows
            for (int row = 0; row < grid.GetLength(1); row++)
            {
                // Loop over the columns
                for (int col = 0; col < grid.GetLength(0); col++)
                {
                    // Only loop over the insides, since the outside nodes are equal to the outside temperature
                    if (col != 0 && col != grid.GetLength(0) - 1 && row != 0 && row != grid.GetLength(1) - 1)
                    {
                        // Assume that the starting temperature is 0
                        int leftTemp = grid[col - 1, row];
                        int topTemp = grid[col, row - 1];
                        int rightTemp = grid[col + 1, row];
                        int bottomTemp = grid[col, row + 1];

                        // Calculate the temperature for the node and store it
                        grid[col, row] = CalcTemp(leftTemp, topTemp, rightTemp, bottomTemp);
                    }
                }
            }
        }

        /*
         * ParTempSolver
         * 
         * Uses a parallel version of the Gauss-Seidell algorithm to find the temperature of each
         * section of a 2D metal plate that is submerged in a solution with a temperature of 0 degrees
         * celsius on each of the sides with the exception of the left which has a temperature of
         * 2000 degrees celsius
         * 
         * Parameters:
         * int[,] grid - The input grid of all the nodes that temperatures will be determined for
         */
        static void ParGaussSeidel(int[,] grid)
        {
            // Keep calculating the temperatures until they stabilize

            // Previous Grid
            // Used to make sure the final temperature has been found
            int[,] prevGrid = new int[grid.GetLength(0), grid.GetLength(1)];

            // Use LINQ to compare the grids (http://stackoverflow.com/questions/12446770/how-to-compare-multidimensional-arrays-in-c-sharp)
            // Only stop when the previous is the same as the current
            while (!(prevGrid.Rank == grid.Rank &&
                     Enumerable.Range(0, prevGrid.Rank).All(dimension => prevGrid.GetLength(dimension) == grid.GetLength(dimension)) &&
                     prevGrid.Cast<int>().SequenceEqual(grid.Cast<int>())))
            {
                // Cannot just assign it straight across because that passes the memory space not the values
                prevGrid = (int[,])grid.Clone();

                // Calculate the temperatures
                ParTempSolver(grid);
            }
        }

        /*
         * ParTempSolver
         * 
         * Finds the temperature for each of the nodes
         * 
         * Parameters:
         * int[,] grid - The input grid of all the nodes that temperatures will be calculated for
         */
        static void ParTempSolver(int[,] grid)
        {
            // Set up the Parallel loop for full parallelism
            ParallelOptions options = new ParallelOptions();
            options.MaxDegreeOfParallelism = -1;

            // Chunk size
            int chunkSize = 8;

            // Divide the array in chunks each rows/cc size.
            Parallel.For(0, chunkSize, options, i => {

                // Number of rows the core will process.
                int lRows = grid.GetLength(1)/chunkSize;

                // Starting row
                int start = i * lRows;

                // Ending row
                int end = (lRows + 1) * i;

                // Prevent IndexOutOfBoundsException for last processor.
                if (end >= grid.GetLength(1))
                {
                    end = grid.GetLength(1);
                }

                // Loop over the rows
                for (int row = start; row < end; row ++ )
                {
                    // Loop over the columns
                    for (int col = 0; col < grid.GetLength(0); col++)
                    {
                        // Only loop over the insides, since the outside nodes are equal to the outside temperature
                        if (col != 0 && col != grid.GetLength(0) - 1 && row != 0 && row != grid.GetLength(1) - 1)
                        {
                            // Assume that the starting temperature is 0
                            int leftTemp = grid[col - 1, row];
                            int topTemp = grid[col, row - 1];
                            int rightTemp = grid[col + 1, row];
                            int bottomTemp = grid[col, row + 1];

                            // Calculate the temperature for the node and store it
                            grid[col, row] = CalcTemp(leftTemp, topTemp, rightTemp, bottomTemp);
                        }
                    }
                }
            });
        }

        /*
         * CalcTemp
         * 
         * Performs the calculation for the node temperature
         * 
         * Parameters:
         * int left   - The left node temperature
         * int top    - The top node temperature
         * int right  - The right node temperature
         * int bottom - The bottom node temperature
         * 
         * Returns: The node temperature
         */
        static int CalcTemp(int left, int top, int right, int bottom)
        {
            // Find the temperature by averaging the four grids around it
            return (left + top + right + bottom) / 4;
        }
    }
}
