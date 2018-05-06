using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace NumericMethods.Blazor.SVD
{
    public static class HivensTransformation
    {
        public static Matrix<double> GenerateRotation(Matrix<double> m, int i, double t)
        {
            var c = 1 / Math.Sqrt(1 + t * t);
            var s = t / Math.Sqrt(1 + t * t);
            var identity = Matrix<double>.Build.DenseIdentity(m.RowCount, m.RowCount);
            identity[i, i] = c;
            identity[i + 1, i + 1] = c;
            identity[i + 1, i] = s;
            identity[i, i + 1] = -s;
            return identity;
        }
    }
}
