using System;
using System.Data;
using System.Linq;


namespace GaussAlgorithm;

public class Solver
{
	public double[] Solve(double[][] matrix, double[] freeMembers)
    {
        int n = matrix.Length; 
        int m = matrix[0].Length;
        double[][] augmentedMatrix = new double[n][];
        for (int i = 0; i < n; i++)
        {
            augmentedMatrix[i] = new double[m + 1];
            for (int j = 0; j < m; j++)
            {
                augmentedMatrix[i][j] = matrix[i][j];
            }
            augmentedMatrix[i][m] = freeMembers[i];
        }
        for(int row = 0; row<n; row++) 
        {
            if (matrix[row].All(x => x <= 1e-10) && freeMembers[row] >= 1e-10)
            {
                throw new NoSolutionException("No solution!");
            }
        }
        var rankMatrix = CalculateRank( matrix);
        var rankAugmentedMatrix = CalculateRank(augmentedMatrix);
        if (rankMatrix < rankAugmentedMatrix)
        {
            throw new NoSolutionException(matrix, freeMembers, augmentedMatrix);
        }
        bool[] usedRows = new bool[n];
        for (int col = 0; col < m; col++)
        {
            int pivotRow = -1;
            for (int row = 0; row < n; row++)
            {
                if (!usedRows[row] && Math.Abs(augmentedMatrix[row][col]) > 1e-10)
                {
                    pivotRow = row;
                    break;
                }
            }

            if (pivotRow == -1)
                continue;

            usedRows[pivotRow] = true;

            double pivotValue = augmentedMatrix[pivotRow][col];
            for (int j = col; j <= m; j++)
            {
                augmentedMatrix[pivotRow][j] /= pivotValue;
            }

            for (int row = 0; row < n; row++)
            {
                if (row != pivotRow)
                {
                    double factor = augmentedMatrix[row][col];
                    for (int j = col; j <= m; j++)
                    {
                        augmentedMatrix[row][j] -= factor * augmentedMatrix[pivotRow][j];
                    }
                }
            }
        }
        double[] solution = new double[m];
        for (int i = 0; i < m; i++)
        {
            solution[i] = 0; 
        }
        for (int row = n - 1; row >= 0; row--)
        {
            int col = -1;
            for (int j = 0; j < m; j++)
            {
                if (Math.Abs(augmentedMatrix[row][j]) > 1e-10)
                {
                    col = j;
                    break;
                }
            }

            if (col != -1)
            {
                solution[col] = augmentedMatrix[row][m];
                for (int j = col + 1; j < m; j++)
                {
                    solution[col] -= augmentedMatrix[row][j] * solution[j];
                }
            }
        }
        return solution;
    }
    static int CalculateRank( double[][] m)
    {
        var matrix = m.Select(row => row.ToArray()).ToArray();
        int rows = matrix.Length;
        int cols = matrix[0].Length;
        int rank = Math.Min(rows, cols); 
        for (int row = 0; row < rank; row++)
        {
            if (matrix[row][row] == 0)
            {
                bool foundNonZero = false;
                for (int i = row + 1; i < rows; i++)
                {
                    if (matrix[i][row] != 0)
                    {
                        SwapRows(matrix, row, i);
                        foundNonZero = true;
                        break;
                    }
                }
                if (!foundNonZero)
                {
                    rank--;
                    if (row < rank)
                    {
                        for (int i = 0; i < rows; i++)
                        {
                            matrix[i][row] = matrix[i][rank];
                        }
                    }
                    row--;
                    continue;
                }
            }

            for (int i = row + 1; i < rows; i++)
            {
                double factor = matrix[i][row] / matrix[row][row];
                for (int j = row; j < cols; j++)
                {
                    matrix[i][j] -= factor * matrix[row][j];
                }
            }
        }

        return rank;
    }

    static void SwapRows(double[][] matrix, int row1, int row2)
    {
        double[] temp = matrix[row1];
        matrix[row1] = matrix[row2];
        matrix[row2] = temp;
    }
}
