using bem.BoundaryValueProblem;

double nu = 0.3;
double kappa = 0.0;

double M_phi(int n, double x, List<double> y) =>
    (1 - nu * nu) / x / x * (n * n * y[0] + x * y[1]) + nu * y[2];

double M_rphi(int n, double x, List<double> y) =>
    -n / x / x * (1 - nu) * (n * n * y[0] + x * y[1]);


double M_rphi_r(int n, double x, List<double> y) =>
    n * x * (-nu * (n * n / x / x * y[0] + y[1] / x) + y[2]) * (1 - nu);

double Q_phi(int n, double x, List<double> y) =>
    1 / x * (n * M_phi(n, x, y) + M_rphi_r(n, x, y));


List<Func<double, List<double>, double>> system(int n)
{
    return new List<Func<double, List<double>, double>> {
    (x, y) => -y[1],
    (x, y) => -nu * (n * n / x / x * y[0] + y[1] / x) + y[2],
    (x, y) => -1 / x * (y[3]-n*M_rphi(n,x,y)-M_phi(n,x,y))+y[3],
    (x, y) => kappa*kappa*y[0]-y[3]+n*Q_phi(n,x,y)// заглушка
    };
};


var _system = system(0);
Dictionary<int, double> left = new Dictionary<int, double> { { 0, 1 }, { 1, 0 } };
Dictionary<int, double> right = new Dictionary<int, double> { { 2, 0 }, { 3, 0 } };
BoundaryValueProblem boundaryValueProblem = new BoundaryValueProblem(_system, left, right, 0.1);
var sol = boundaryValueProblem.Solve();