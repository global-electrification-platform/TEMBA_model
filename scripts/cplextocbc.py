"""Converts a CPLEX solution file into the same format produced by CBC

"""
import sys
from typing import List, Union


class ConvertLine(object):
    """Abstract class which defines the interface to the family of convertors

    Inherit this class and implement the ``_do_it()`` method to produce the
    data to be written out into a new format

    Example
    -------
    >>> cplex_line = "AnnualCost	REGION	CDBACKSTOP	1.0	0.0	137958.8400384134"
    >>> convertor = RegionTechnology()
    >>> convertor.convert()
    VariableName(REGION,TECHCODE01,2015)       42.69         0\n
    VariableName(REGION,TECHCODE01,2017)       137958.84         0\n
    """

    def __init__(self, data: List):
        self.data = data

    def _do_it(self) -> List:
        raise NotImplementedError()

    def convert(self) -> List[str]:
        """Perform the conversion
        """
        cbc_data = []
        variable, dimensions, values = self._do_it()

        for index, value in enumerate(values):

            year = 2015 + index
            if (value not in ["0.0", "0", ""]) and (year <= 2070):

                try:
                    value = float(value)
                except ValueError:
                    value = 0

                full_dims = ",".join(dimensions + (str(year),))

                formatted_data = "0 {0}({1}) {2} 0\n".format(
                    variable,
                    full_dims,
                    value
                    )

                cbc_data.append(formatted_data)

        return cbc_data


class RegionTimeSliceTechnologyMode(ConvertLine):

    def _do_it(self) -> List:
        """Produces output indexed by Region, Timeslice, Tech and Mode

        ``VariableName(REGION,SD1D,TECHCODE01,2,2015)       42.69         0\n``

        """
        variable = self.data[0]
        region = self.data[1]
        timeslice = self.data[2]
        technology = self.data[3]
        mode = self.data[4]
        values = self.data[5:]

        dimensions = (region, timeslice, technology, mode)

        return (variable, dimensions, values)


class RegionTechnology(ConvertLine):

    def _do_it(self) -> List:
        """Produces output indexed by dimensions Region and Technology

        ``0\tVariableName(REGION,TECHCODE01,2015)\t42.69\t0\n``

        """
        variable = self.data[0]
        region = self.data[1]
        technology = self.data[2]

        dimensions = (region, technology)

        values = self.data[3:]

        return (variable, dimensions, values)


def process_line(line: str) -> Union[List, None]:
    """Processes an individual line in a CPLEX file

    A different ConvertLine implementation is chosen depending upon the
    variable name

    Arguments
    ---------
    line : str
    """
    row_as_list = line.split('\t')
    variable = row_as_list[0]
    if variable in ['NewCapacity',
                    'TotalCapacityAnnual',
                    'CapitalInvestment',
                    'AnnualFixedOperatingCost',
                    'AnnualVariableOperatingCost']:
        convertor = RegionTechnology(row_as_list)
    elif variable in ['RateOfActivity']:
        convertor = RegionTimeSliceTechnologyMode(row_as_list)
    else:
        convertor = None

    return convertor


def convert_cplex_file(cplex_filename, output_filename):
    """Converts a CPLEX solution file into that of the CBC solution file

    Arguments
    ---------
    cplex_filename : str
        Path to the transformed CPLEX solution file
    output_filename : str
        Path for the processed data to be written to
    """

    with open(output_filename, 'w') as cbc_file:
        with open(cplex_filename, 'r') as cplex_file:
            for linenum, line in enumerate(cplex_file):
                convertor = process_line(line)
                try:
                    if convertor:
                        data = convertor.convert()
                        cbc_file.writelines(data)
                except ValueError:
                    msg = "Error caused at line {}: {}"
                    raise ValueError(msg.format(linenum, line))


class TestCplexRead:

    def test_different_format(self):

        fixture = "AnnualFixedOperatingCost	REGION	AOBACKSTOP	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0		"

        expected = []

        convertor = process_line(fixture)
        actual = convertor.convert()
        assert actual == expected

    def test_read_in_line(self):

        fixture = "AnnualFixedOperatingCost	REGION	CDBACKSTOP	0.0	0.0	137958.8400384134	305945.38410619126	626159.9611543404	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0"

        expected = [
                    "AnnualFixedOperatingCost(REGION,CDBACKSTOP,2017)                             137958.8                      0",
                    "AnnualFixedOperatingCost(REGION,CDBACKSTOP,2018)                            305945.38                      0",
                    "AnnualFixedOperatingCost(REGION,CDBACKSTOP,2019)                            626159.96                      0"]

        convertor = process_line(fixture)
        actual = convertor.convert()
        assert actual == expected


    def test_rate_of_activity(self):

        fixture = """RateOfActivity	REGION	S1D1	CGLFRCFURX	1	0.0	0.0	0.0	0.0	0.0	0.3284446367303371	0.3451714779880536	0.3366163200621617	0.3394945166233896	0.3137488154250392	0.28605725055560716	0.2572505015401749	0.06757558148965725	0.0558936625751148	0.04330608461292407	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0"""

        expected = [
            "RateOfActivity(REGION,S1D1,CGLFRCFURX,1,2020)                              0.32844464                      0",
            "RateOfActivity(REGION,S1D1,CGLFRCFURX,1,2021)                              0.34517148                      0",
            "RateOfActivity(REGION,S1D1,CGLFRCFURX,1,2022)                              0.33661632                      0",
            "RateOfActivity(REGION,S1D1,CGLFRCFURX,1,2023)                              0.33949452                      0",
            "RateOfActivity(REGION,S1D1,CGLFRCFURX,1,2024)                              0.31374882                      0",
            "RateOfActivity(REGION,S1D1,CGLFRCFURX,1,2025)                              0.28605725                      0",
            "RateOfActivity(REGION,S1D1,CGLFRCFURX,1,2026)                               0.2572505                      0",
            "RateOfActivity(REGION,S1D1,CGLFRCFURX,1,2027)                             0.067575581                      0",
            "RateOfActivity(REGION,S1D1,CGLFRCFURX,1,2028)                             0.055893663                      0",
            "RateOfActivity(REGION,S1D1,CGLFRCFURX,1,2029)                             0.043306085                      0"]

        convertor = process_line(fixture)
        actual = convertor.convert()
        assert actual == expected


if __name__ == '__main__':

    if len(sys.argv) != 3:
        msg = "Usage:\npython {} <cplex_file> <output_file>\n"
        print(msg.format(sys.argv[0]))
        sys.exit(1)

    cplex_file = sys.argv[1]
    cbc_file = sys.argv[2]

    convert_cplex_file(cplex_file, cbc_file)