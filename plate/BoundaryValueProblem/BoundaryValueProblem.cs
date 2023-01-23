using RungeKutta;

namespace bem.BoundaryValueProblem
{
    public class BoundaryValueProblem
    {
        // однородные уравнения системы
        private readonly List<Func<double, List<double>, double>> _equations;

        // неоднородная часть системы
        private readonly Dictionary<int, Func<double, double>> _inhomogeneous;

        // краевые условия на левом крае
        private readonly Dictionary<int, double> _leftConditions;

        // условия на правом крае
        private readonly Dictionary<int, double> _rightConditions;

        // левый конец отрезка
        private readonly double _left_boundary;

        // точность расчёта
        private readonly double _epsilon;

        /// <summary>
        /// метод создаёт набор начальных условия для вспомогательных задач Коши
        /// </summary>
        /// <returns>
        /// набор начальных условия для вспомогательных задач Коши в виде 
        /// матрицы вещественных чисел 
        /// </returns>
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

        /// <summary>
        /// решения вспомогательных задач Коши на правом конце отрезка
        /// </summary>
        /// <returns>
        /// матрица вещественных чисел
        /// </returns>
        private List<List<double>> CauchyProblemSolutions()
        {
            var conditions = InitialConditions();
            var size = conditions.Count;
            List<List<double>> solutions = new List<List<double>>();
            OdeSolver odeSolver = new OdeSolver(InhomogeneousSystem(), ButcherTableau.Rkf78, _epsilon);
            var solution = odeSolver.Solve(_left_boundary, 1.0, conditions.First());
            solutions.Add(solution);
            odeSolver = new OdeSolver(_equations, ButcherTableau.Rkf78, _epsilon);
            for (int i = 1; i < size; i++)
            {
                List<double> initials = new List<double>(_equations.Count);
                for (int ii = 0; ii < _equations.Count; ii++)
                {
                    initials.Add(conditions[i][ii]);
                }
                solution = odeSolver.Solve(_left_boundary, 1.0, initials);
                solutions.Add(solution);
            }
            return solutions;
        }

        /// <summary>
        /// построение вектора начальных условий для решения краевой задачи
        /// </summary>
        /// <param name="solutions">
        /// решения вспомогательных задач Коши на правом конце отрезка
        /// </param>
        /// <returns></returns>
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

        /// <summary>
        /// решение краевой задачи
        /// </summary>
        /// <returns>
        /// решение краевой задачи на правой конце отрезка
        /// </returns>
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

        /// <summary>
        /// система для решения неоднородной задачи
        /// </summary>
        /// <returns>система для решения неоднородной задачи</returns>
        private List<Func<double, List<double>, double>> InhomogeneousSystem()
        {
            List<Func<double, List<double>, double>> inhomogeneousSystem =
                new List<Func<double, List<double>, double>>(_equations);
            foreach (var item in _inhomogeneous)
            {
                inhomogeneousSystem[item.Key]
                    = (x, y) => _equations[item.Key](x, y) + item.Value(x);
            }
            return inhomogeneousSystem;
        }

        /// <summary>
        /// конструктор
        /// </summary>
        /// <param name="equations">уравнения</param>
        /// <param name="leftConditions">краевые условия на правом конце отрезка</param>
        /// <param name="rightConditions">краевые условия на левом концен отрезка</param>
        /// <param name="inner">левая граница отрезка</param>
        /// <param name="epsilon">точность вычислений</param>
        public BoundaryValueProblem(List<Func<double, List<double>, double>> equations,
            Dictionary<int, double> leftConditions, Dictionary<int, double> rightConditions,
            Dictionary<int, Func<double, double>> inhomogeneous,
            double inner = 0, double epsilon = 0.1e-6)
        {
            _equations = equations;
            _leftConditions = leftConditions;
            _rightConditions = rightConditions;
            _inhomogeneous = inhomogeneous;
            _left_boundary = inner;
            _epsilon = epsilon;
        }
    }
}
