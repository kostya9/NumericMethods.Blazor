using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericMethods.Blazor.SvdDecomposition.Householder
{
    /// <summary>
    /// Performs householder transformations to change the matrix to a up two-diagonal one
    /// </summary>
    public class HigherHouseholderTransformation : IHouseholderTransformation
    {
        private readonly Matrix<double> _matrix;
        private List<Matrix<double>> _rightHouseHolder;
        private List<Matrix<double>> _leftHouseHolder;
        private List<Matrix<double>> _hivens;


        public Matrix<double> Result { get; private set; }
        public IReadOnlyCollection<Matrix<double>> RightHouseholder => _rightHouseHolder;
        public IReadOnlyCollection<Matrix<double>> LeftHouseHolder => _leftHouseHolder;
        public IReadOnlyCollection<Matrix<double>> Hivens => _hivens;


        public HigherHouseholderTransformation(Matrix<double> matrix)
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

            if (rows == cols)
            {
                leftHCount = t - 1;
                rightHCount = t - 2;
            }
            else if (rows > cols)
            {
                leftHCount = t;
                rightHCount = t - 2;
            }
            else
            {
                leftHCount = t - 1;
                rightHCount = t - 1;

                if (cols > rows + 1)
                    rightHCount = t;
            }

            Matrix<double> matrix = _matrix;
            int leftIteration = 0;
            int rightIteration = 0;
            bool performed = true;
            while (performed)
            {
                performed = false;
                if (leftIteration < leftHCount)
                {
                    var leftMatrix = GenerateLeftMatrix(matrix, leftIteration++);
                    matrix = leftMatrix * matrix;
                    _leftHouseHolder.Add(leftMatrix);
                    performed = true;
                }

                if (rightIteration < rightHCount)
                {
                    var rightMatrix = GenerateRightMatrix(matrix, rightIteration++);
                    matrix = matrix * rightMatrix;
                    _rightHouseHolder.Add(rightMatrix);
                    performed = true;
                }
            }

            _hivens = new List<Matrix<double>>();
            if (cols > rows)
            {
                var colToRace = rows;
                for (int r = rows - 1; r >= 0; r--)
                {
                    var tan = -matrix[r, colToRace] / matrix[r, r];
                    var hivens = HivensTransformation.GenerateRotation(matrix, r, rows, tan, cols);
                    _hivens.Add(hivens);
                    matrix = matrix * hivens;
                }
            }

            Result = matrix;
        }

        private Matrix<double> GenerateLeftMatrix(Matrix<double> matrix, int iteration)
        {
            var i = iteration;
            var element = matrix[i, i];
            var sign = Math.Sign(-element);
            var beta = sign * Math.Sqrt(Enumerable.Range(i, matrix.RowCount - i).Select(r => matrix[r, i] * matrix[r, i]).Sum());
            var mu = 1 / Math.Sqrt(2 * beta * beta - 2 * beta * element);
            var w = mu * Vector<double>.Build.DenseOfEnumerable(Enumerable.Range(0, matrix.RowCount).Select(r =>
            {
                if (r < i)
                    return 0;
                if (r == i)
                    return element - beta;
                return matrix[r, i];
            }));

            return Matrix<double>.Build.DenseIdentity(matrix.RowCount, matrix.RowCount) - 2 * w.ToColumnMatrix() * w.ToRowMatrix();
        }

        private Matrix<double> GenerateRightMatrix(Matrix<double> matrix, int iteration)
        {
            var i = iteration;
            var element = matrix[i, i + 1];
            var sign = Math.Sign(-element);
            var beta = sign * Math.Sqrt(Enumerable.Range(i + 1, matrix.ColumnCount - i - 1).Select(c => matrix[i, c] * matrix[i, c]).Sum());
            var mu = 1 / Math.Sqrt(2 * beta * beta - 2 * beta * element);
            var w = mu * Vector<double>.Build.DenseOfEnumerable(Enumerable.Range(0, matrix.ColumnCount).Select(c =>
            {
                if (c < i + 1)
                    return 0;
                if (c == i + 1)
                    return element - beta;
                return matrix[i, c];
            }));

            return Matrix<double>.Build.DenseIdentity(matrix.ColumnCount, matrix.ColumnCount) - 2 * w.ToColumnMatrix() * w.ToRowMatrix();
        }
    }
}
