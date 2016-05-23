using System;
using System.Diagnostics; // For the stopwatch
using System.Threading.Tasks; // For the Parallel
using System.Linq; // For the 2D array comparison
using System.Drawing; // For the color
using System.Drawing.Imaging; // For the image output

namespace FalseColoring
{
    class Program
    {

        static void Main(string[] args)
        {
            // Set up the grid dimensions
            // Parallel will be slower for smaller grids
            int width = 500;
            int length = 500;

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

            // Output the results
            //PrintGrid(grid);
            OutputImage(grid, "image1.bmp"); // Image stored in the bin/Debug folder of the project directory

            // Fill the 2D array with the starting values
            SetupGrid(grid);

            // Solve in parallel
            sw.Reset();
            sw.Start();
            ParGaussSeidel(grid);
            sw.Stop();
            Console.WriteLine("Elapsed time taken parallel: {0}", sw.ElapsedMilliseconds);

            // Output the results
            //PrintGrid(grid);
            OutputImage(grid, "image2.bmp"); // Image stored in the bin/Debug folder of the project directory

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

        /*
         * AssignColor
         * 
         * Assigns a color to the given node temperature
         * 
         * Parameters:
         * int temp - The node temperature
         * 
         * Returns: The corresponding color
         */
        static Color AssignColor(int temp)
        {
            double[] colorVals = new double[3];

            // Color map is as follows:
            // 0    -> (0, 0, 128)
            // 500  -> (0, 255, 0)
            // 1000 -> (255, 255, 0)
            // 1500 -> (255, 128, 0)
            // 2000 -> (255, 0, 0)
            // Source: http://dsp.stackexchange.com/a/4679

            // If the temp is not on the map, use linear interpolation to 
            // determine the color for the temperature
            if (temp == 0)
            {
                colorVals[0] = 0;
                colorVals[1] = 0;
                colorVals[2] = 128;
            }
            else if (temp == 500)
            {
                colorVals[0] = 0;
                colorVals[1] = 255;
                colorVals[2] = 0;
            }
            else if (temp == 1000)
            {
                colorVals[0] = 255;
                colorVals[1] = 255;
                colorVals[2] = 0;
            }
            else if (temp == 1500)
            {
                colorVals[0] = 255;
                colorVals[1] = 128;
                colorVals[2] = 0;
            }
            else if (temp == 2000)
            {
                colorVals[0] = 255;
                colorVals[1] = 0;
                colorVals[2] = 0;
            }
            else
            {
                int[] leftRGB = new int[3];
                int[] rightRGB = new int[3];
                double multiplier = 0.0;

                // Determine the left and right RGB values for the calculation
                if (temp > 0 && temp < 500)
                {
                    leftRGB[0] = 0;
                    leftRGB[1] = 0;
                    leftRGB[2] = 128;
                    rightRGB[0] = 0;
                    rightRGB[1] = 255;
                    rightRGB[2] = 0;
                    multiplier = 0.25;
                }
                else if (temp > 500 && temp < 1000)
                {
                    leftRGB[0] = 0;
                    leftRGB[1] = 255;
                    leftRGB[2] = 0;
                    rightRGB[0] = 255;
                    rightRGB[1] = 255;
                    rightRGB[2] = 0;
                    multiplier = 0.5;
                }
                else if (temp > 1000 && temp < 1500)
                {
                    leftRGB[0] = 255;
                    leftRGB[1] = 255;
                    leftRGB[2] = 0;
                    rightRGB[0] = 255;
                    rightRGB[1] = 128;
                    rightRGB[2] = 0;
                    multiplier = 0.75;
                }
                else if (temp > 1500 && temp < 2000)
                {
                    leftRGB[0] = 255;
                    leftRGB[1] = 128;
                    leftRGB[2] = 0;
                    rightRGB[0] = 255;
                    rightRGB[1] = 0;
                    rightRGB[2] = 0;
                    multiplier = 1.0;
                }

                // Assign the colors
                for (int i = 0; i < 3; i++ )
                {
                    colorVals[i] = ((double) temp / 2000.0) * (double) leftRGB[i] / multiplier + ((2000.0 - (double) temp) / 2000.0) * rightRGB[i];

                    // Make sure the colors don't exceed the accepted ranges
                    if (colorVals[i] > 255.0) colorVals[i] = 255.0;
                }
            }

            return Color.FromArgb((int) Math.Ceiling(colorVals[0]), (int) Math.Ceiling(colorVals[1]), (int) Math.Ceiling(colorVals[2]));
        }

        static void PrintGrid(int[,] grid)
        {
            // Print the grid
            for (int row = 0; row < grid.GetLength(1); row++)
            {
                for (int col = 0; col < grid.GetLength(0); col++)
                {
                    if (col == grid.GetLength(0) - 1)
                    {
                        System.Console.WriteLine(grid[col, row]);
                    }
                    else
                    {
                        System.Console.Write(grid[col, row] + "-");
                    }
                }
            }
        }

        static void OutputImage(int[,] grid, string imageName)
        {
            // Set up the image
            Bitmap bitmap = new Bitmap(grid.GetLength(0), grid.GetLength(1));

            // Add the pixels to the image
            for (int row = 0; row < grid.GetLength(1); row++)
            {
                for (int col = 0; col < grid.GetLength(0); col++)
                {
                    bitmap.SetPixel(col, row, AssignColor(grid[col, row]));
                }
            }

            // Save it
            bitmap.Save(imageName);
        }
    }
}
