using MathNet.Numerics.LinearAlgebra;

namespace NumericMethods.Blazor.SvdDecomposition.SvdDecomposition
{
    public interface ISvdDecomposition
    {
        Matrix<double> U { get; }
        Matrix<double> Sigma { get; }
        Matrix<double> V { get; }

        void Perform();
    }
}
