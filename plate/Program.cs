using bem;
using bem.BoundaryValueProblem;
using RungeKutta;


class Plate
{
    static double nu = 0.3;
    //static double kappa = 20.0;
    static double innerRadius = 0.00001;

    static List<Func<double, List<double>, double>> system(int n, double kappa)
    {
        return new List<Func<double, List<double>, double>> {
    (x, y) => -y[1],
    (x, y) => -nu*y[1]/x+y[2],
    (x, y) => (1-nu*nu)*y[1]/x/x-(1-nu)*y[2]/x+y[3],
    (x, y) => -y[3]/x - kappa*kappa*y[0]
    };
    }

    /// <summary>
    /// амплитудно-частотная характеристика
    /// </summary>
    /// <param name="maxKappa">наибольшая частота (построение происходит с нуля)</param>
    /// <param name="numPoints">количество узловых точек</param>
    /// <param name="func">функция (АЧХ можно задать аналитическим выражением)</param>
    /// <returns>АЧХ в виде словаря</returns>
    static Dictionary<double, double> FrequencyResponse(double maxKappa,
        int numPoints, Func<double, double> func)
    {
        Dictionary<double, double> result = new Dictionary<double, double>();
        var h = maxKappa / numPoints;
        for (int i = 0; i < numPoints; i++)
        {
            var kappa = i * h;
            result.Add(kappa, func(kappa));
        }
        return result;
    }

    /// <summary>
    /// перегрузка FrequenceResponse с минимальной частотой
    /// </summary>
    /// <param name="minKappa"></param>
    /// <param name="maxKappa"></param>
    /// <param name="numPoints"></param>
    /// <param name="func"></param>
    /// <returns></returns>    
    static Dictionary<double, double> FrequencyResponse(double minKappa,
        double maxKappa, int numPoints, Func<double, double> func)
    {
        Dictionary<double, double> result = new Dictionary<double, double>();
        var h = (maxKappa - minKappa) / numPoints;
        for (int i = 0; i < numPoints; i++)
        {
            var kappa = minKappa + i * h;
            result.Add(kappa, func(kappa));
        }
        return result;
    }

    /// <summary>
    /// частотное уравнение
    /// </summary>
    /// <param name="kappa">максимальная частота</param>
    /// <returns>ачх в виде списка вещественных чисел</returns>
    static double FrecquencyEquation(double kappa)
    {
        var equations = system(0, kappa);
        OdeSolver odeSolver = new OdeSolver(equations, ButcherTableau.RungeKuttaFeldberg);
        var solution1 = odeSolver.Solve(innerRadius, 1.0, new List<double> { 1, 0, 0, 0 });
        var solution2 = odeSolver.Solve(innerRadius, 1.0, new List<double> { 0, 0, 0, 1 });
        Console.WriteLine(kappa);
        return solution1[2] * solution2[3] - solution1[3] * solution2[2];
    }


    /// <summary>
    /// частотное уравнение
    /// </summary>
    /// <param name="kappa">максимальная частота</param>
    /// <returns>ачх в виде списка вещественных чисел</returns>
    static double FrecquencyResponse(double kappa)
    {
        var equations = system(0, kappa);
        Dictionary<int, double> left = new Dictionary<int, double> { { 0, 1 }, { 1, 0 } };
        Dictionary<int, double> right = new Dictionary<int, double> { { 2, 0 }, { 3, 0 } };
        BoundaryValueProblem boundaryValueProblem = new BoundaryValueProblem(equations, left, right, new Dictionary<int, Func<double, double>>(), innerRadius);
        Console.WriteLine(kappa);
        return boundaryValueProblem.Solve()[0];
    }


    //OdeSolver odeSolver = new OdeSolver(system(0), ButcherTableau.RungeKuttaFeldberg);
    //var solution1 = odeSolver.Solve(innerRadius, 1.0, new List<double> { 0, 0, 1, 0 });
    static void Main()
    {
        var fr = FrequencyResponse(50.0, 50, x => FrecquencyResponse(x));
        Tikz.Plot(fr, "plot.txt");
    }
}