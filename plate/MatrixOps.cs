using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace bem
{
    static class MatrixOps
    {
        /// <summary>
        /// Решение системы линейных алгебраических уравнений
        /// </summary>
        /// <param name="matrix">матрица</param>
        /// <param name="rightPart">правая часть</param>
        /// <returns></returns>
        public static double[] Solve(double[,] matrix, double[] rightPart)
        {
            if (matrix.GetLength(0) != rightPart.Length)
            {
                throw new ArgumentException("");
            }
            var lu = LUDecompose(matrix, out var p);
            var rp = TransformTheRightPart(rightPart, p);
            return ThereAndBackAgain(lu, rp);
        }

        /// <summary>
        /// Решение СЛАУ методом Гаусса
        /// </summary>
        /// <param name="lu">LU-разложение исходной матрицы</param>
        /// <param name="rp">правая часть</param>
        /// <returns>решение</returns>
        private static double[] ThereAndBackAgain(double[,] lu, double[] rp)
        {
            double[] y = new double[rp.Length];
            for (int i = 0; i < rp.Length; i++)
            {
                double sum = 0;
                for (int ii = 0; ii < i; ii++)
                {
                    sum += y[ii] * lu[i, ii];
                }
                y[i] = rp[i] - sum;
            }

            for (int i = rp.Length - 1; i >= 0; i--)
            {
                double sum = 0;
                for (int ii = i + 1; ii < rp.Length; ii++)
                {
                    sum += y[ii] * lu[i, ii];
                }
                y[i] = (y[i] - sum) / lu[i, i];
            }
            return y;
        }

        /// <summary>
        /// Преобразование правой части СЛАУ (при произведении LU-разложения
        /// строки матрицы меняются местами, правую часть следует изменить
        /// соответственно)
        /// </summary>
        /// <param name="rightPart">првая часть</param>
        /// <param name="p">массив хранит  информацию о перестановках строк</param>
        /// <returns>новая правая часть</returns>
        static double[] TransformTheRightPart(double[] rightPart, int[] p)
        {
            if (rightPart.Length != p.Length)
            {
                throw new ArgumentException("");
            }
            double[] rp = new double[rightPart.Length];
            for (int i = 0; i < rp.Length; i++)
            {
                rp[i] = rightPart[p[i]];
            }
            return rp;
        }

        /// <summary>
        /// LU-разложение 
        /// </summary>
        /// <param name="matrix">матрица СЛАУ</param>
        /// <param name="p">массив с информацией о перестановках</param>
        /// <returns>матрица LU</returns>
        static double[,] LUDecompose(double[,] matrix, out int[] p)
        {
            if (matrix.GetLength(0) != matrix.GetLength(1))
            {
                throw new ArgumentException("");
            }
            var size = matrix.GetLength(0);
            var elems = (double[,])matrix.Clone();
            p = Enumerable.Range(0, size).ToArray();
            for (int i = 0; i < size; i++)
            {
                var imax = Imax(i, elems, size);
                if (imax != i)
                {
                    SwapRows(p, elems, i, imax);
                }
                for (int ii = i + 1; ii < size; ii++)
                {
                    elems[ii, i] /= elems[i, i];
                    for (int k = i + 1; k < size; k++)
                    {
                        elems[ii, k] -= elems[ii, i] * elems[i, k];
                    }
                }
            }
            return elems;
        }

        /// <summary>
        /// одновременная перестановка строк в матрице и в векторе p
        /// </summary>
        /// <param name="p">вектор перестановок</param>
        /// <param name="elems">копия матрицы СЛАУ</param>
        /// <param name="i">номер строки</param>
        /// <param name="j">номер строки</param>
        private static void SwapRows(int[] p, double[,] elems, int i, int j)
        {
            (p[i], p[j]) = (p[j], p[i]);
            for (int ii = 0; ii < elems.GetLength(1); ii++)
            {
                (elems[i, ii], elems[j, ii]) = (elems[j, ii], elems[i, ii]);
            }
        }

        private static int Imax(int i, double[,] elems, int size)
        {
            var imax = i;
            var max = Math.Abs(elems[i, i]);
            for (int ii = i; ii < size; ii++)
            {
                var abs = Math.Abs(elems[ii, i]);
                if (abs > max)
                {
                    max = abs;
                    imax = ii;
                }
            }
            if (max < double.Epsilon)
            {
                throw new ArgumentException("матрица вырождена");
            }

            return imax;//(imax, max);
        }

        public static double[,] CommonMatrix(double[,] stiffnessMatrix, double[,] massMatrix, double kappa)
        {
            var newMass = Multiply(kappa * kappa, massMatrix);
            return SumMatrix(stiffnessMatrix, newMass);
        }

        public static double[,] SumMatrix(double[,] stiffnessMatrix, double[,] newMass)
        {
            var result = (double[,])stiffnessMatrix.Clone();
            for (int i = 0; i < result.GetLength(0); i++)
            {
                for (int ii = 0; ii < result.GetLength(1); ii++)
                {
                    result[i, ii] -= newMass[i, ii];
                }
            }
            return result;
        }

        public static double[,] Multiply(double kappa, double[,] massMatrix)
        {
            var result = (double[,])massMatrix.Clone();
            for (int i = 0; i < result.GetLength(0); i++)
            {
                for (int ii = 0; ii < result.GetLength(1); ii++)
                {
                    result[i, ii] *= kappa;
                }
            }
            return result;
        }

        public static double[] Multiply(double[,] matrix, double[] vector)
        {
            if (matrix.GetLength(1) != vector.Length)
            {
                throw new ArgumentOutOfRangeException("");
            }
            double[] result = new double[matrix.GetLength(0)];
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int ii = 0; ii < matrix.GetLength(1); ii++)
                {
                    result[i] += matrix[i, ii] * vector[ii];
                }
            }
            return result;
        }

        public static double[] MultiplyT(double[,] matrix, double[] vector)
        {
            if (matrix.GetLength(0) != vector.Length)
            {
                throw new ArgumentOutOfRangeException("");
            }
            double[] result = new double[matrix.GetLength(1)];
            for (int i = 0; i < matrix.GetLength(1); i++)
            {
                for (int ii = 0; ii < matrix.GetLength(0); ii++)
                {
                    result[i] += matrix[ii, i] * vector[ii];
                }
            }
            return result;
        }

        public static double[] SumVectors(double[] left, double[] right)
        {
            if (left.Length != right.Length)
            {
                throw new ArgumentOutOfRangeException("");
            }
            double[] result = new double[left.Length];
            for (int i = 0; i < left.Length; i++)
            {
                result[i] = left[i] + right[i];
            }
            return result;
        }

        /// <summary>
        /// норма
        /// </summary>
        /// <param name="vs"></param>
        /// <returns></returns>
        public static double Norm(this IEnumerable<double> vs)
        {
            return Math.Sqrt(vs.Aggregate(0.0, (x, y) => x + y * y));
        }

        /// <summary>
        /// нормализация
        /// </summary>
        /// <param name="vs"></param>
        /// <returns></returns>
        public static double[] Normalize(this IEnumerable<double> vs)
        {
            var enumerable = vs as double[] ?? vs.ToArray();
            var norm = enumerable.Norm();
            if (Math.Abs(norm) < double.Epsilon)
            {
                return new[] { 1.0 };
            }
            return enumerable.Select(x => x / norm).ToArray();
        }

        public static double[] Normalize(this IEnumerable<double> vs, double norm)
        {
            var enumerable = vs as double[] ?? vs.ToArray();
            if (Math.Abs(norm) < double.Epsilon)
            {
                return new[] { 1.0 };
            }
            return enumerable.Select(x => x / norm).ToArray();
        }
    }
}
