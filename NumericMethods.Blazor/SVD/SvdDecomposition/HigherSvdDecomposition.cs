using MathNet.Numerics.LinearAlgebra;
using NumericMethods.Blazor.SvdDecomposition.SvdDecomposition;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace NumericMethods.Blazor.SVD.SvdDecomposition
{
    public class HigherSvdDecomposition : ISvdDecomposition
    {

        public HigherSvdDecomposition(Matrix<double> m, Matrix<double> l, Matrix<double> r)
        {

        }

        public Matrix<double> U => Matrix<double>.Build.DiagonalIdentity(5);

        public Matrix<double> Sigma => Matrix<double>.Build.DiagonalIdentity(5);

        public Matrix<double> V => Matrix<double>.Build.DiagonalIdentity(5);

        public void Perform()
        {
            // TODO: Implement
        }
    }
}
