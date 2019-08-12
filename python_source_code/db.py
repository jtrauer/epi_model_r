
from sqlalchemy import create_engine
import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
from python_source_code.curve import scale_up_function


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
        self.tabs_of_interest = ["BCG", "Aggregated estimates"]

    def update_csv_reads(self, input_path="../xls/*.csv"):
        """
        load csvs from input_path
        """
        csv_file_list = glob.glob(input_path)
        for filename in csv_file_list:
            data_frame = pd.read_csv(filename)
            data_frame.to_sql(filename.split("\\")[1].split(".")[0], con=self.engine, if_exists="replace")

    def update_xl_reads(self, input_path="../xls/*.xlsx"):
        """
        load excel spreadsheet from input_path
        """
        excel_file_list = glob.glob(input_path)
        for filename in excel_file_list:
            xls = pd.ExcelFile(filename)

            # for single tab in spreadsheet
            if len(xls.sheet_names) == 1:
                df_name = xls.sheet_names[0]
                df = pd.read_excel(filename, sheet_name=df_name)
                df.to_sql(df_name, con=self.engine, if_exists="replace")
                self.output_to_user("now reading '%s' tab of '%s' file" % (df_name, filename))

            # if multiple tabs
            else:
                for n_sheets, sheet in enumerate(xls.sheet_names):
                    if sheet in self.tabs_of_interest:
                        header_3_sheets = ["rate_birth_2015", "life_expectancy_2015"]
                        n_header = 3 if sheet in header_3_sheets else 0
                        df = pd.read_excel(filename, sheet_name=sheet, header=n_header)
                        self.output_to_user("now reading '%s' tab of '%s' file" % (sheet, filename))

                        # to read constants and time variants
                        if sheet == "constants":
                            sheet = filename.replace(".xlsx", "").split("_")[1] + "_constants"
                        if sheet == "time_variants":
                            sheet = filename.replace(".xlsx", "").split("_")[1] + "_time_variants"
                        df.to_sql(sheet, con=self.engine, if_exists="replace")

    def output_to_user(self, comment):
        """
        report progress to user if requested
        """
        if self.verbose:
            print(comment)

    def db_query(self, table_name, is_filter="", value="", column="*"):
        """
        method to query table_name
        """
        query = "Select %s from  %s" % (column, table_name)
        if is_filter and value:
            query = query + " Where %s = \'%s\'" % (is_filter, value)
        return pd.read_sql_query(query, con=self.engine)


if __name__ == "__main__":

    # standard code to update the database
    input_database = InputDB(verbose=True)
    input_database.update_xl_reads()
    input_database.update_csv_reads()

    # example of accessing once loaded
    res = input_database.db_query("gtb_2016", column="c_cdr", is_filter="country", value="Mongolia")
    cdr_mongolia = res["c_cdr"].values
    res = input_database.db_query("gtb_2016", column="year", is_filter="country", value="Mongolia")
    cdr_mongolia_year = res["year"].values
    spl = scale_up_function(cdr_mongolia_year, cdr_mongolia, smoothness=0.2, method=5)
    times = list(np.linspace(1950, 2014, 1e3))
    scaled_up_cdr = []
    for t in times:
        scaled_up_cdr.append(spl(t))
    plt.plot(cdr_mongolia_year, cdr_mongolia, "ro", times, scaled_up_cdr)
    plt.title("CDR from GTB 2015")
    plt.show()
