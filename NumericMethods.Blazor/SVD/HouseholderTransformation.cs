using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace NumericMethods.Blazor.SVD
{
    /// <summary>
    /// Performs householder transformations to change the matrix to a down two-diagonal one
    /// </summary>
    public class HouseholderTransformation
    {
        private readonly Matrix<double> _matrix;
        private List<Matrix<double>> _rightHouseHolder;
        private List<Matrix<double>> _leftHouseHolder;

        public Matrix<double> Result { get; private set; }
        public IReadOnlyCollection<Matrix<double>> RightHouseholder => _rightHouseHolder;
        public IReadOnlyCollection<Matrix<double>> LeftHouseHolder => _leftHouseHolder;

        public HouseholderTransformation(Matrix<double> matrix)
        {
            _matrix = matrix;
        }

        public void Perform()
        {
            var cols = _matrix.ColumnCount;
            var rows = _matrix.RowCount;
            var t = Math.Min(rows, cols);

            _leftHouseHolder = new List<Matrix<double>>();
            int leftHCount = 0;
            _rightHouseHolder = new List<Matrix<double>>();
            int rightHCount = 0;

            if(rows == cols)
            {
                leftHCount = t - 2;
                rightHCount = t - 1;
            }
            else if(rows > cols)
            {
                leftHCount = t - 1;
                rightHCount = t - 1;

                if (rows > cols + 1)
                    leftHCount = t;
            }
            else
            {
                leftHCount = t - 2;
                rightHCount = t;
            }


            Matrix<double> matrix = _matrix;
            int leftIteration = 0;
            int rightIteration = 0;
            bool performed = true;
            while(performed)
            {
                performed = false;
                if (rightIteration < rightHCount)
                {
                    var rightMatrix = GenerateRightMatrix(matrix, rightIteration++);
                    matrix = matrix * rightMatrix;
                    _rightHouseHolder.Add(rightMatrix);
                    performed = true;
                }

                if (leftIteration < leftHCount)
                {
                    var leftMatrix = GenerateLeftMatrix(matrix, leftIteration++);
                    matrix = leftMatrix * matrix;
                    _leftHouseHolder.Add(leftMatrix);
                    performed = true;
                }
            }

            if(rows > cols)
            {
                var rowToRace = cols;
                for(int c = cols - 1; c >= 0; c--)
                {
                    var tan = - matrix[rowToRace, c] / matrix[c, c];
                    var hivens = HivensTransformation.GenerateRotation(matrix, cols, c, tan);
                    matrix = hivens * matrix;
                }
            }

            Result = matrix;
        }

        private Matrix<double> GenerateLeftMatrix(Matrix<double> matrix, int iteration)
        {
            var i = iteration;
            var element = matrix[i + 1, i];
            var sign = Math.Sign(-element);
            var beta = sign * Math.Sqrt(Enumerable.Range(i + 1, matrix.RowCount - i - 1).Select(r => matrix[r, i] * matrix[r, i]).Sum());
            var mu = 1 / Math.Sqrt(2 * beta * beta - 2 * beta * element);
            var w = mu * Vector<double>.Build.DenseOfEnumerable(Enumerable.Range(0, matrix.RowCount).Select(r =>
            {
                if (r < i + 1)
                    return 0;
                if (r == i + 1)
                    return element - beta;
                return matrix[r, i];
            }));

            return Matrix<double>.Build.DenseIdentity(matrix.RowCount, matrix.RowCount) - 2 * w.ToColumnMatrix() * w.ToRowMatrix();
        }

        private Matrix<double> GenerateRightMatrix(Matrix<double> matrix, int iteration)
        {
            var i = iteration;
            var element = matrix[i, i];
            var sign = Math.Sign(-element);
            var beta = sign * Math.Sqrt(Enumerable.Range(i, matrix.ColumnCount - i).Select(r => matrix[i, r] * matrix[i, r]).Sum());
            var mu = 1 / Math.Sqrt(2 * beta * beta - 2 * beta * element);
            var w = mu * Vector<double>.Build.DenseOfEnumerable(Enumerable.Range(0, matrix.ColumnCount).Select(c =>
            {
                if (c < i)
                    return 0;
                if (c == i)
                    return element - beta;
                return matrix[i, c];
            }));

            return Matrix<double>.Build.DenseIdentity(matrix.ColumnCount, matrix.ColumnCount) - 2 * w.ToColumnMatrix() * w.ToRowMatrix();
        }
    }
}
