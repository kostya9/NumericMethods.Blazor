using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace NumericMethods.Blazor.SvdDecomposition.Householder
{
    public interface IHouseholderTransformation
    {
        Matrix<double> Result { get; }
        IReadOnlyCollection<Matrix<double>> RightHouseholder { get; }
        IReadOnlyCollection<Matrix<double>> LeftHouseHolder { get; }
        IReadOnlyCollection<Matrix<double>> Hivens { get; }

        void Perform();
    }
}
