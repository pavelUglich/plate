namespace bem
{
    class Tikz
    {

        public static void FileHeader(string fileName, out StreamWriter streamWriter)
        {
            streamWriter = new StreamWriter(fileName);
            streamWriter.Write("\\begin{tikzpicture}[scale=2.0]\\tiny\n");
            streamWriter.Write("\\begin{axis}[grid]\n");
        }

        public static void FileFooter(StreamWriter streamWriter)
        {
            streamWriter.Write("\\end{axis}\n");
            streamWriter.Write("\\end{tikzpicture}\n");
            streamWriter.Close();
        }

        public static void Plot(Dictionary<string, IEnumerable<double>> curves, string fileName)
        {
            FileHeader(fileName, out var streamWriter);
            foreach (var curve in curves)
            {
                AddTheCurve(curve.Value, streamWriter, curve.Key);
            }
            FileFooter(streamWriter);
        }


        /// <summary>
        /// график на отрезке от 0 до 1
        /// </summary>
        /// <param name="numbers"></param>
        public static void Plot(IEnumerable<double> numbers, string fileName)
        {
            FileHeader(fileName, out var streamWriter);
            AddTheCurve(numbers, streamWriter);
            FileFooter(streamWriter);
        }

        public static void AddTheCurve(IEnumerable<double> numbers,
            StreamWriter streamWriter, string color = "black", double a = 0, double b = 1)
        {
            streamWriter.Write("\\addplot[line width = 0.25mm, smooth, ");
            streamWriter.Write($"{color}] plot coordinates{{\n");
            var array = numbers.ToArray();
            var h = (b - a) / (array.Length - 1);
            for (int i = 0; i < array.Length; i++)
            {
                var x = (a + i * h).ToString().Replace(',', '.');
                var y = array[i].ToString().Replace(',', '.');
                streamWriter.Write($"({x}, {y})");
            }
            streamWriter.Write("};\n");
        }


        public static void AddTheCurve(Dictionary<double, double> numbers,
            StreamWriter streamWriter, string color = "black")
        {
            streamWriter.Write("\\addplot[line width = 0.25mm, smooth, ");
            streamWriter.Write($"{color}] plot coordinates{{\n");
            foreach (var item in numbers)
            {
                var x = item.Key.ToString().Replace(',', '.');
                var y = item.Value.ToString().Replace(',', '.');
                streamWriter.Write($"({x}, {y})");
            }
            streamWriter.Write("};\n");
        }


        /// <summary>
        /// график на отрезке от 0 до 1
        /// </summary>
        /// <param name="numbers"></param>
        public static void Plot(Dictionary<double, double> numbers, string fileName)
        {
            FileHeader(fileName, out var streamWriter);
            streamWriter.Write("\\addplot[line width = 0.25mm, smooth");
            streamWriter.Write("] plot coordinates{\n");
            foreach (var item in numbers)
            {
                var x = item.Key.ToString().Replace(',', '.');
                var y = item.Value.ToString().Replace(',', '.');
                streamWriter.Write($"({x}, {y})");
            }
            streamWriter.Write("};\n");
            FileFooter(streamWriter);
        }


        public static void Plot(Dictionary<double, List<double>> numbers,
            int component, string fileName)
        {
            FileHeader(fileName, out var streamWriter);
            streamWriter.Write("\\addplot[line width = 0.25mm, smooth");
            streamWriter.Write("] plot coordinates{\n");
            foreach (var item in numbers)
            {
                var x = item.Key.ToString().Replace(',', '.');
                var y = item.Value[component].ToString().Replace(',', '.');
                streamWriter.Write($"({x}, {y})");
            }
            streamWriter.Write("};\n");
            FileFooter(streamWriter);
        }


        /// <summary>
        /// график на отрезке от 0 до 1
        /// </summary>
        /// <param name="numbers"></param>
        public static void Plot(Dictionary<double, double> numbers, StreamWriter streamWriter)
        {
            streamWriter.Write("\\addplot[line width = 0.25mm, smooth");
            streamWriter.Write("] plot coordinates{\n");
            foreach (var item in numbers)
            {
                var x = item.Key.ToString().Replace(',', '.');
                var y = item.Value.ToString().Replace(',', '.');
                streamWriter.Write($"({x}, {y})");
            }
            streamWriter.Write("};\n");
        }
    }
}
