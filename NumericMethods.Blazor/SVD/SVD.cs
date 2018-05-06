using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace NumericMethods.Blazor.SVD
{
    public class SVD
    {
        private readonly Matrix<double> _matrix;
        private readonly Matrix<double> _leftH;
        private readonly Matrix<double> _rightH;

        public Matrix<double> U { get; private set; }
        public Matrix<double> Sigma { get; private set; }
        public Matrix<double> V { get; private set; }

        public SVD(Matrix<double> matrix, Matrix<double> leftH, Matrix<double> rightH)
        {
            if (matrix.ColumnCount != matrix.RowCount)
                throw new ArgumentException("Expected matrix to be quadratic");

            var epsilon = Math.Pow(10, -4);
            for (int r = 0; r < matrix.RowCount; r++)
            {
                for (int c = 0; c < matrix.ColumnCount; c++)
                {
                    if (r != c && r != c + 1)
                        if (Math.Abs(matrix[r, c]) > epsilon)
                            throw new ArgumentException($"Expected a twodiagonal matrix. Error at m[{r}, {c}]={matrix[r, c]}");
                }
            }

            _matrix = matrix;
            _leftH = leftH;
            _rightH = rightH;
        }

        public void Perform()
        {
            var t = GenerateFirstT();
            Console.WriteLine(t);
            Sigma = t * _matrix;
        }

        private Matrix<double> GenerateFirstT()
        {
            var n = _matrix.RowCount - 1;
            var d = (_matrix[n, n] * _matrix[n, n] - _matrix[n - 1, n - 1] * _matrix[n - 1, n - 1] + _matrix[n - 1, n] * _matrix[n - 1, n] - _matrix[n - 2, n - 1]) / 
                (2 * _matrix[n - 1, n] * _matrix[n - 1, n - 1]);
            var r = d >= 0 ? (-d - Math.Sqrt(1 + d * d)) : (d + Math.Sqrt(1 + d * d));
            var thau = _matrix[n, n] * _matrix[n, n] + _matrix[n - 1, n] - _matrix[n - 1, n] * _matrix[n - 1, n - 1] / r;
            var t = _matrix[0, 0] * _matrix[0, 1] / (_matrix[0, 0] * _matrix[0, 0] - thau);
            Console.WriteLine(t);
            return GenerateHivens(_matrix, 0, t);
        }

        private Matrix<double> GenerateT(Matrix<double> b, int iteration)
        {
            var i = iteration;
            var t = b[i - 1, i + 1] / b[i, i];
            return GenerateHivens(b, iteration, t);
        }

        private Matrix<double> GenerateS(Matrix<double> p, int iteration)
        {
            var i = iteration;
            var t = p[i + i, i] / p[i, i];
            return GenerateHivens(p, iteration, t);
        }

        private Matrix<double> GenerateHivens(Matrix<double> m, int i, double t)
        {
            var c = 1 / Math.Sqrt(1 + t * t);
            var s = t / Math.Sqrt(1 + t * t);
            var identity = Matrix<double>.Build.DiagonalIdentity(m.RowCount);
            identity[i, i] = c;
            identity[i + 1, i + 1] = c;
            identity[i + 1, i] = s;
            identity[i, i + 1] = -s;
            return identity;
        }
    }
}
