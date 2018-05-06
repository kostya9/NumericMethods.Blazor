using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace NumericMethods.Blazor.Import
{
    public class CsvImport
    {
        public static event EventHandler<CsvImport> OnImport;

        public int[,] Result { get; }

        public CsvImport(int[,] importData)
        {
            Result = importData;
        }

        public static void FromFile(string base64)
        {
            byte[] bytes = Convert.FromBase64String(base64);

            using (var reader = new MemoryStream(bytes))
            {
                FromFile(reader);
            }
        }

        private static void FromFile(Stream stream)
        {
            List<int[]> values = new List<int[]>();
            using (TextReader streamReader = new StreamReader(stream))
            {
                using (var reader = new CsvHelper.CsvParser(streamReader))
                {
                    for(var record = reader.Read(); record != null; record = reader.Read())
                    {
                        values.Add(record.Select(int.Parse).ToArray());
                    }

                    int[,] matrix = new int[values.Count(), values.First().Count()];

                    for (var r = 0; r < matrix.GetLength(0); r++)
                    {
                        for (var c = 0; c < matrix.GetLength(1); c++)
                        {
                            matrix[r, c] = values[r][c];
                        }
                    }

                    OnImport?.Invoke(null, new CsvImport(matrix));
                }
            }
        }
    }
}
