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

            var n = _matrix.RowCount - 1;
            var firstT = GenerateFirstT();
            var b  =  _matrix * firstT;

            int tCount = n - 1;
            int sCount = n - 1;
            int tIteration = 1;
            int sIteration = 0;
            bool performed = true;
            while(performed)
            {
                performed = false;
                if(sIteration < sCount)
                {
                    var s = GenerateS(b, sIteration++);
                    b = s * b;
                    performed = true;
                }

                if(tIteration < tCount)
                {
                    var t = GenerateT(b, tIteration++);
                    b =  b * t;
                    performed = true;
                }
                Console.WriteLine(b);
            }
            Sigma = b;
        }

        private Matrix<double> GenerateFirstT()
        {
            var n = _matrix.RowCount - 1;
            var d = (Square(_matrix[n, n]) - Square(_matrix[n - 1, n - 1]) + Square(_matrix[n, n - 1]) - Square(_matrix[n - 1, n - 2])) / 
                (2 * _matrix[n, n - 1] * _matrix[n - 1, n - 1]);
            var r = d >= 0 ? (-d - Math.Sqrt(1 + d * d)) : (-d + Math.Sqrt(1 + d * d));
            var thau = Square(_matrix[n, n]) + Square(_matrix[n, n - 1]) - (_matrix[n, n - 1] * _matrix[n - 1, n - 1]) / r;
            var t = _matrix[0, 0] * _matrix[1, 0] / (Square(_matrix[0, 0]) - thau);
            return HivensTransformation.GenerateRotation(_matrix, 0, t);
        }

        private Matrix<double> GenerateT(Matrix<double> b, int iteration)
        {
            var i = iteration;
            var t = b[i - 1, i + 1] / b[i, i];
            return HivensTransformation.GenerateRotation(b, iteration, t);
        }

        private Matrix<double> GenerateS(Matrix<double> p, int iteration)
        {
            var i = iteration;
            var t = p[i, i + 1] / p[i, i];
            return HivensTransformation.GenerateRotation(p, iteration, t);
        }

        private double Square(double d)
        {
            return d * d;
        }
    }
}
