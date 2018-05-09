using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace NumericMethods.Blazor.SVD
{
    /// <summary>
    /// SVD for lower diagonal matrix
    /// </summary>
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

        // TODO: Doesn't work as expected
        public void Perform()
        {
            List<Matrix<double>> sMatrices = new List<Matrix<double>>();
            List<Matrix<double>> tMatrices = new List<Matrix<double>>();

            Sigma = Iterate(_matrix);
        }

        private Matrix<double> Iterate(Matrix<double> m)
        {
            var n = m.RowCount;
            var firstS = GenerateFirstS();
            Console.WriteLine("firstS=" + firstS);
            var b = firstS * m;

            Console.WriteLine(b);

            int tCount = n - 1;
            int sCount = n - 1;
            int tIteration = 0;
            int sIteration = 1;
            bool performed = true;
            while (performed)
            {
                performed = false;
                if (tIteration < tCount)
                {
                    var t = GenerateT(b, tIteration++);
                    b = b * t;
                    performed = true;
                    Console.WriteLine(b);
                }

                if (sIteration < sCount)
                {
                    var s = GenerateS(b, sIteration++);
                    b = s * b;
                    performed = true;
                    Console.WriteLine(b);
                }
            }

            return b;
        }

        private Matrix<double> GenerateFirstS()
        {
            var n = _matrix.RowCount - 1;
            var d = (Square(_matrix[n, n]) - Square(_matrix[n - 1, n - 1]) + Square(_matrix[n - 1, n]) - Square(_matrix[n - 2, n - 1])) / 
                (2 * _matrix[n - 1, n] * _matrix[n - 1, n - 1]);
            var r = d >= 0 ? (-d - Math.Sqrt(1 + d * d)) : (-d + Math.Sqrt(1 + d * d));
            var thau = Square(_matrix[n, n]) + Square(_matrix[n - 1, n]) - (_matrix[n - 1, n] * _matrix[n - 1, n - 1]) / r;
            var t = _matrix[1, 0] * _matrix[0, 0] / (Square(_matrix[0, 0]) - thau);
            Console.WriteLine("St=" + t);
            return HivensTransformation.GenerateRotation(_matrix, 0, t);
        }

        private Matrix<double> GenerateT(Matrix<double> b, int iteration)
        {
            var i = iteration;
            var t = b[i, i + 1] / b[i, i];
            Console.WriteLine("Tt=" + t); 
            return HivensTransformation.GenerateRotation(b, i, t);
        }

        private Matrix<double> GenerateS(Matrix<double> p, int iteration)
        {
            var i = iteration;
            var t = p[i - 1, i + 1] / p[i - 1, i];
            Console.WriteLine("St=" + t);
            return HivensTransformation.GenerateRotation(p, iteration, t);
        }

        private double Square(double d)
        {
            return d * d;
        }
    }
}
