
from sqlalchemy import create_engine
import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
from python_source_code.curve import scale_up_function


def get_bcg_coverage(database, country_iso_code):
    """
    extract bcg coverage from inputs database

    :param database: sql database
        the database containing the bcg data
    :param country_iso_code: string
        three letter ISO3 code for the country of interest
    :return: bcg_coverage
        pandas data frame with columns years and one row containing the values of BCG coverage in that year
    """
    _bcg_coverage = database.db_query("BCG", is_filter="ISO_code", value=country_iso_code)
    _bcg_coverage = _bcg_coverage.filter(items=[column for column in _bcg_coverage.columns if column.isdigit()])
    return {int(key): value / 1e2 for key, value in zip(list(_bcg_coverage.columns), _bcg_coverage.loc[0, :])
            if value is not None}


def get_all_iso3_from_bcg(database):
    """
    check which iso3 country codes are available from the bcg database

    :param database: sql database
        the database containing the bcg data
    :return: list
        all the iso3 strings available from the bcg database
    """

    return database.db_query("bcg", column="ISO_code")["ISO_code"].tolist()


class InputDB:
    """
    methods for loading input xls files
    """

    def __init__(self, database_name="../databases/inputs.db", verbose=False):
        """
        initialise sqlite database
        """
        self.database_name = database_name
        self.engine = create_engine("sqlite:///" + database_name, echo=False)
        self.verbose = verbose
        self.headers_lookup = \
            {"../xls/WPP2019_FERT_F03_CRUDE_BIRTH_RATE.xlsx": 16,
             "../xls/WPP2019_F01_LOCATIONS.xlsx": 16,
             "../xls/WPP2019_MORT_F04_1_DEATHS_BY_AGE_BOTH_SEXES.xlsx": 16,
             "../xls/WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx": 16,
             "../xls/life_expectancy_2015.xlsx": 3,
             "../xls/rate_birth_2015.xlsx": 3}
        self.tab_of_interest = \
            {"../xls/WPP2019_FERT_F03_CRUDE_BIRTH_RATE.xlsx": "ESTIMATES",
             "../xls/WPP2019_MORT_F04_1_DEATHS_BY_AGE_BOTH_SEXES.xlsx": "ESTIMATES",
             "../xls/WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx": "ESTIMATES",
             "../xls/WPP2019_F01_LOCATIONS.xlsx": "Location",
             "../xls/coverage_estimates_series.xlsx": "BCG",
             "../xls/gtb_2015.xlsx": "gtb_2015",
             "../xls/gtb_2016.xlsx": "gtb_2016",
             "../xls/life_expectancy_2015.xlsx": "life_expectancy_2015",
             "../xls/rate_birth_2015.xlsx": "rate_birth_2015"}
        self.output_name = \
            {"../xls/WPP2019_FERT_F03_CRUDE_BIRTH_RATE.xlsx": "crude_birth_rate",
             "../xls/WPP2019_MORT_F04_1_DEATHS_BY_AGE_BOTH_SEXES.xlsx": "absolute_deaths",
             "../xls/WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx": "total_population",
             "../xls/WPP2019_F01_LOCATIONS.xlsx": "un_iso3_map",
             "../xls/coverage_estimates_series.xlsx": "bcg",
             "../xls/gtb_2015.xlsx": "gtb_2015",
             "../xls/gtb_2016.xlsx": "gtb_2016",
             "../xls/life_expectancy_2015.xlsx": "life_expectancy_2015",
             "../xls/rate_birth_2015.xlsx": "rate_birth_2015"}
        self.map_df = None

    def update_csv_reads(self, input_path="../xls/*.csv"):
        """
        load csvs from input_path
        """
        csv_file_list = glob.glob(input_path)
        for filename in csv_file_list:
            data_frame = pd.read_csv(filename)
            data_frame.to_sql(filename.split("\\")[1].split(".")[0], con=self.engine, if_exists="replace")

    def update_xl_reads(self, sheets_to_read=glob.glob("../xls/*.xlsx")):
        """
        load excel spreadsheet from input_path

        :param sheets_to_read: iterable
            paths of the spreadsheets to read, which have to be strictly coded in the format suggested above
        """
        for available_file in sheets_to_read:
            filename = "../xls/" + available_file[7: -5] + ".xlsx"
            header_row = self.headers_lookup[filename] if filename in self.headers_lookup else 0
            data_title = self.output_name[filename] if filename in self.output_name else filename
            current_data_frame = pd.read_excel(
                pd.ExcelFile(filename), header=header_row, index_col=1, sheet_name=self.tab_of_interest[filename])
            self.output_to_user("now reading '%s' tab of '%s' file" % (self.tab_of_interest[filename], filename))
            current_data_frame.to_sql(data_title, con=self.engine, if_exists="replace")

    def output_to_user(self, comment):
        """
        report progress to user if requested

        :param comment: str
            string to be output to the user
        """
        if self.verbose:
            print(comment)

    def db_query(self, table_name, is_filter="", value="", column="*"):
        """
        method to query table_name
        """
        query = "Select %s from %s" % (column, table_name)
        if is_filter and value:
            query = query + " Where %s = \'%s\'" % (is_filter, value)
        return pd.read_sql_query(query, con=self.engine)

    def add_iso_to_table(self, table_name):
        """
        add the mapped iso3 code to a table that only contains the un country code

        :param table_name: str
            name of the spreadsheet to perform this on
        """
        self.get_un_iso_map()
        table_with_iso = pd.merge(
            self.db_query(table_name=table_name), self.map_df, left_on='Country code', right_on='Location code')
        if "Index" in table_with_iso.columns:
            table_with_iso = table_with_iso.drop(columns=["Index"])
        table_with_iso.to_sql(table_name + "_mapped", con=self.engine, if_exists="replace")

    def get_un_iso_map(self):
        """
        create dictionary structure to map from un three numeric digit codes to iso3 three alphabetical digit codes
        """
        self.map_df = self.db_query(table_name='un_iso3_map')[['Location code', 'ISO3 Alpha-code']].dropna()


if __name__ == "__main__":

    # standard code to update the database
    input_database = InputDB()
    # input_database.update_xl_reads()
    # input_database.update_csv_reads()
    # input_database.add_iso_to_table("crude_birth_rate")
    # input_database.add_iso_to_table("absolute_deaths")
    # input_database.add_iso_to_table("total_population")

    # example of accessing once loaded
    # times = list(np.linspace(1950, 2020, 1e3))
    # extract data for BCG vaccination for a particular country
    # for country in get_all_iso3_from_bcg(input_database):
    #     bcg_coverage = get_bcg_coverage(input_database, country)
    #     if len(bcg_coverage) == 0:
    #         print("no BCG vaccination data available for %s" % country)
    #         continue
    #     print("plotting BCG vaccination data and fitted curve for %s" % country)
    #     bcg_coverage_function = scale_up_function(
    #           bcg_coverage.keys(), bcg_coverage.values(), smoothness=0.2, method=5)
    #     plt.plot(list(bcg_coverage.keys()), list(bcg_coverage.values()), "ro")
    #     plt.plot(times, [bcg_coverage_function(time) for time in times])
    #     plt.title(country)
    #     plt.show()
