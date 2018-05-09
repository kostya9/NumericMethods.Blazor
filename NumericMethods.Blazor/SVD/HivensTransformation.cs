using MathNet.Numerics.LinearAlgebra;
using System;

namespace NumericMethods.Blazor.SvdDecomposition
{
    public static class HivensTransformation
    {
        public static Matrix<double> GenerateRotation(Matrix<double> m, int i, int k, double tan, int size)
        {
            var c = 1 / Math.Sqrt(1 + tan * tan);
            var s = c * tan;
            var identity = Matrix<double>.Build.DenseIdentity(size, size);
            identity[i, i] = c;
            identity[k, k] = c;
            identity[i, k] = s;
            identity[k, i] = -s;
            return identity;
        }

        public static Matrix<double> GenerateRotation(Matrix<double> m, int i, double t)
            => GenerateRotation(m, i, i + 1, t, m.RowCount);
    }
}
