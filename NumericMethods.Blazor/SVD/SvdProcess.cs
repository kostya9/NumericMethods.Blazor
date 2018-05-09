using MathNet.Numerics.LinearAlgebra;
using NumericMethods.Blazor.SVD.SvdDecomposition;
using NumericMethods.Blazor.SvdDecomposition.Householder;
using NumericMethods.Blazor.SvdDecomposition.SvdDecomposition;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace NumericMethods.Blazor.SVD
{
    public class SvdProcess
    {
        private readonly Func<Matrix<double>, IHouseholderTransformation> _householderFactory;
        private readonly Func<Matrix<double>, Matrix<double>, Matrix<double>, ISvdDecomposition> _decompositionFactory;

        private SvdProcess(Func<Matrix<double>, IHouseholderTransformation> householderTransformation, Func<Matrix<double>, Matrix<double>, Matrix<double>, ISvdDecomposition> decomposition)
        {
            _householderFactory = householderTransformation;
            _decompositionFactory = decomposition;
        }

        public SvdProcessOutput Perform(Matrix<double> matrix)
        {
            var householder = _householderFactory(matrix);
            householder.Perform();
            var left = householder.LeftHouseHolder.Reverse().Aggregate(Matrix<double>.Build.DiagonalIdentity(matrix.RowCount), (a, c) => a * c);
            var right = householder.RightHouseholder.Aggregate(Matrix<double>.Build.DiagonalIdentity(matrix.ColumnCount), (a, c) => a * c);
            var result = householder.Result;
            if (result.RowCount > result.ColumnCount)
            {
                var times = result.RowCount - result.ColumnCount;
                while (times-- != 0)
                    result = result.RemoveRow(result.RowCount - 1);

            }
            else if (result.ColumnCount > result.RowCount)
            {
                var times = result.ColumnCount - result.RowCount;
                while (times-- != 0)
                    result = result.RemoveColumn(result.ColumnCount - 1);

            }
            var svd = _decompositionFactory(result, left, right);
            svd.Perform();
            return new SvdProcessOutput
            {
                LeftHouseholders = householder.LeftHouseHolder,
                RightHouseholders = householder.RightHouseholder,
                RaceHivens = householder.Hivens,
                HouseholderResult = householder.Result,
                Iterations = new List<SvdProcessOutput.SvdIteration>()
                {
                    new SvdProcessOutput.SvdIteration()
                    {
                        U = svd.U,
                        Sigma = svd.Sigma,
                        V = svd.V
                    }
                }
            };
        }

        public static SvdProcess CreateLowerTwoDiagonal()
        {
            return new SvdProcess((m) => new LowerHouseholderTransformation(m), (m, l, r) => new LowerSvdDecomposition(m, l, r));
        }

        public static SvdProcess CreateHigherTwoDiagonal()
        {
            return new SvdProcess((m) => new HigherHouseholderTransformation(m), (m, l, r) => new HigherSvdDecomposition(m, l, r));
        }
    }

    public class SvdProcessOutput
    {
        public IReadOnlyCollection<Matrix<double>> LeftHouseholders { get; set; }
        public IReadOnlyCollection<Matrix<double>> RightHouseholders { get; set; }
        public IReadOnlyCollection<Matrix<double>> RaceHivens { get; set; }
        public Matrix<double> HouseholderResult { get; set; }
        public IReadOnlyCollection<SvdIteration> Iterations { get; set; }

        public class SvdIteration
        {
            public IReadOnlyCollection<Matrix<double>> S { get; set; }
            public IReadOnlyCollection<Matrix<double>> T { get; set; }
            public Matrix<double> U { get; set; }
            public Matrix<double> Sigma { get; set; }
            public Matrix<double> V { get; set; }
        }
    }
}
