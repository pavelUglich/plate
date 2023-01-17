using System;
using System.Collections.Generic;
using System.Linq;
using RungeKutta;

namespace bem.BoundaryValueProblem
{
    public class BoundaryValueProblem
    {
        private readonly List<Func<double, List<double>, double>> _equations;
        private readonly Dictionary<int, double> _leftConditions;
        private readonly Dictionary<int, double> _rightConditions;
        private readonly double _epsilon;


        private List<List<double>> InitialConditions()
        {
            var size = _equations.Count;
            List<List<double>> result = new List<List<double>>(size + 1);
            List<double> first = new List<double>(_equations.Count);
            for (int i = 0; i < _equations.Count; i++)
            {
                var newValue = _leftConditions.ContainsKey(i) ? _leftConditions[i] : 0;
                first.Add(newValue);
            }
            result.Add(first);
            for (int i = 0; i < size; i++)
            {
                if (_leftConditions.ContainsKey(i)) continue;
                var zeros = Enumerable.Repeat(0.0, _equations.Count).ToList();
                zeros[i] = 1.0;
                result.Add(zeros);
            }
            return result;
        }

        private List<List<double>> CauchyProblemSolutions()
        {
            var conditions = InitialConditions();
            var size = conditions.Count;
            List<List<double>> solutions = new List<List<double>>();
            OdeSolver odeSolver = new OdeSolver(_equations, ButcherTableau.Rkf78, _epsilon);
            for (int i = 0; i < size; i++)
            {
                List<double> initials = new List<double>(_equations.Count);
                for (int ii = 0; ii < _equations.Count; ii++)
                {
                    initials.Add(conditions[i][ii]);
                }
                var solution = odeSolver.Solve(0.0, 1.0, initials);
                solutions.Add(solution);
            }
            return solutions;
        }

        private double[] GetTheInitialConditions(IReadOnlyList<List<double>> solutions)
        {
            var size = _rightConditions.Count;
            double[,] matrix = new double[size, size];
            double[] rightPart = new double[size];
            int counter = 0;
            foreach (var condition in _rightConditions)
            {
                rightPart[counter] = condition.Value - solutions[0][condition.Key];
                for (int i = 0; i < size; i++)
                {
                    matrix[counter, i] = solutions[i + 1][condition.Key];
                }
                ++counter;
            }
            return MatrixOps.Solve(matrix, rightPart);
        }

        public List<double> Solve()
        {
            var solutions = CauchyProblemSolutions();
            var initialConditions = GetTheInitialConditions(solutions);
            List<double> result = solutions[0];
            for (int i = 1; i < solutions.Count; i++)
            {
                for (int ii = 0; ii < solutions[i].Count; ii++)
                {
                    result[ii] += initialConditions[i - 1] * solutions[i][ii];
                }
            }
            return result;
        }

        public BoundaryValueProblem(List<Func<double, List<double>, double>> equations,
            Dictionary<int, double> leftConditions, Dictionary<int, double> rightConditions, double epsilon)
        {
            _equations = equations;
            _leftConditions = leftConditions;
            _rightConditions = rightConditions;
            _epsilon = epsilon;
        }
    }
}
