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
        // Set the temperature of the leftmost nodes
        private static int maxTemp = 10000;

        static void Main(string[] args)
        {
            // Set up the grid dimensions
            int width = 100;
            int length = 200;

            // Declare the 2D arrays
            int[,] gridSeq = new int[length, width];
            int[,] gridPar = new int[length, width];

            // Fill the grids with the starting values
            SetupGrid(gridSeq);
            SetupGrid(gridPar);

            // Keep track of the run time
            Stopwatch sw = new Stopwatch();

            // Solve sequentially
            Console.WriteLine("Computing Sequentially . . .");
            sw.Start();
            GaussSeidel(gridSeq, false);
            sw.Stop();
            //PrintGrid(gridSeq);
            Console.WriteLine("Elapsed Time (Sequential): {0}\n", sw.ElapsedMilliseconds);
            OutputImage(gridSeq, "imageSeq.bmp"); // Image stored in the bin/Debug folder of the project directory

            // Solve in parallel
            Console.WriteLine("Computing in Parallel . . .");
            sw.Reset();
            sw.Start();
            GaussSeidel(gridPar, true);
            sw.Stop();
            //PrintGrid(gridPar);
            Console.WriteLine("Elapsed Time (Parallel): {0}\n", sw.ElapsedMilliseconds);
            OutputImage(gridPar, "imagePar.bmp"); // Image stored in the bin/Debug folder of the project directory

            // Program is done, make sure to the let user see the results
            OpenImage("imagePar.bmp");
            OpenImage("imageSeq.bmp");
            System.Console.WriteLine("Press Any Key to Exit");
            Console.ReadKey();

        }

        /*
         * SetupGrid
         * 
         * Sets the starting temperature for each of the nodes assuming that the leftmost column
         * has a temperature of maxTemp degrees celsius and everything else has a temperature of
         * 0 degrees celsius
         * 
         * Parameters:
         * int[,] grid - The input grid of all the nodes
         * int maxTemp - The temperature of the leftmost column
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
                        // Have to initialize the left side to the maximum temperature
                        grid[0, row] = maxTemp;
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
         * GaussSeidel
         * 
         * Uses the Gauss-Seidell algorithm to find the temperatures of each section of a 2D 
         * metal plate that is submerged in a solution with a temperature of 0 degrees
         * celsius on each of the sides with the exception of the left which has a temperature of
         * maxTemp degrees celsius
         * 
         * Parameters:
         * int[,] grid - The input grid of all the nodes that temperatures will be determined for
         */
        static void GaussSeidel(int[,] grid, bool parallel)
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
                if (parallel == true) ParTempSolver(grid);
                else SeqTempSolver(grid);
            }
        }

        /*
         * SeqTempSolver
         * 
         * Finds the temperature for each of the nodes sequentially
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
         * Finds the temperature for each of the nodes in parallel
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
            int cc = Environment.ProcessorCount;

            // Divide the array in chunks each rows/cc size.
            Parallel.For(0, cc, options, i => {

                // Number of rows the core will process.
                int lRows = grid.GetLength(1) / cc;

                // Starting row
                int start = i * lRows;

                // Ending row
                int end = lRows * (i+1);

                // Prevent IndexOutOfBoundsException for last processor.
                if(i == cc - 1)
                {
                    end = grid.GetLength(1);
                }

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
            // Stores the RGB values
            double[] colorVals = new double[3];

            // Begin building the color map
            // Based off of information from http://dsp.stackexchange.com/a/4679
            // and https://en.wikipedia.org/wiki/Linear_interpolation

            // Color markers (cold = blue, lukewarm = white, hot = red)
            if (temp == 0)
            {
                colorVals[0] = 0;
                colorVals[1] = 0;
                colorVals[2] = 255;
            }
            else if (temp == maxTemp / 2)
            {
                colorVals[0] = 255;
                colorVals[1] = 255;
                colorVals[2] = 255;
            }
            else if (temp == maxTemp)
            {
                colorVals[0] = 255;
                colorVals[1] = 0;
                colorVals[2] = 0;
            }
            else
            {
                // If the temp is not on the map, use linear interpolation to 
                // determine the color for the temperature

                // Store the RGB values for the color markers the temp is between
                int[] leftRGB = new int[3];
                int[] rightRGB = new int[3];

                // Store the temperature boundary points
                double leftTempRange = 0;
                double rightTempRange = 0;

                // Determine where the temperature fits between the color markers
                if (temp < maxTemp / 2)
                {
                    // Between blue and white
                    leftRGB[0] = 0;
                    leftRGB[1] = 0;
                    leftRGB[2] = 255;
                    rightRGB[0] = 255;
                    rightRGB[1] = 255;
                    rightRGB[2] = 255;
                    leftTempRange = 0.0;
                    rightTempRange = (double) (maxTemp / 2);
                }
                else
                {
                    // Between the white and red
                    leftRGB[0] = 255;
                    leftRGB[1] = 255;
                    leftRGB[2] = 255;
                    rightRGB[0] = 255;
                    rightRGB[1] = 51;
                    rightRGB[2] = 0;
                    leftTempRange = (double)(maxTemp / 2);
                    rightTempRange = (double) maxTemp;
                }

                // Assign the color based on the temperature using linear interpolation
                for (int i = 0; i < 3; i++ )
                {
                    colorVals[i] = (double) leftRGB[i] + ((double) rightRGB[i] - (double) leftRGB[i]) * (((double) temp - leftTempRange) / (rightTempRange - leftTempRange));
                }
            }

            // Return the RGB values as a Color
            return Color.FromArgb((int) Math.Ceiling(colorVals[0]), (int) Math.Ceiling(colorVals[1]), (int) Math.Ceiling(colorVals[2]));
        }

        /*
         * PrintGrid
         * 
         * Prints the grid to the console. Does not work well for grids that are 
         * larger than the console window.
         * 
         * Parameters:
         * int[,] grid - input grid to be printed
         */
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

        /*
         * OutputImage
         * 
         * Outputs the grid as a BMP false color image.
         * 
         * Parameters:
         * int[,] grid - The input grid to be turned into a false color image
         * string imageName - The filename of the outputted image
         */
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

            // Save the image
            bitmap.Save(imageName);
        }

        /*  
         * OpenImage
         *  
         * Opens the image for viewing. Image will be in running directory of program.
         *  
         * Parameters:
         * string filename - The name of the file that will be opened
         */
        static void OpenImage(string filename)
        {
            // Based off of information from http://stackoverflow.com/a/1283593
            Process process = new Process();
            process.EnableRaisingEvents = false;
            process.StartInfo.FileName = filename;
            process.Start();
        }
    }
}
