using MathNet.Numerics.LinearAlgebra;
using System;

namespace NumericMethods.Blazor.SVD
{
    public static class HivensTransformation
    {
        public static Matrix<double> GenerateRotation(Matrix<double> m, int i, int k, double t)
        {
            var c = 1 / Math.Sqrt(1 + t * t);
            var s = c * t;
            var identity = Matrix<double>.Build.DenseIdentity(m.RowCount, m.RowCount);
            identity[i, i] = c;
            identity[k, k] = c;
            identity[i, k] = s;
            identity[k, i] = -s;
            return identity;
        }

        public static Matrix<double> GenerateRotation(Matrix<double> m, int i, double t)
            => GenerateRotation(m, i, i + 1, t);
    }
}
