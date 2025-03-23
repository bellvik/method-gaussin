using System;
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

        bool[] usedRows = new bool[n];
        for (int col = 0; col < m; col++)
        {
            int pivotRow = -1;
            for (int row = 0; row < n; row++)
            {
                if (!usedRows[row] && Math.Abs(augmentedMatrix[row][col]) > 0)
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
                if (Math.Abs(augmentedMatrix[row][j]) > 0)
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
}