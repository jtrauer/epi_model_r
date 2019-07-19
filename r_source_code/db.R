library(RSQLite)
library(XLConnect)


load_csv = function(input_path="..//xls", database_name="..//databases//Inputs.db"){
    # load csvs from input_path
    conn <- dbConnect(SQLite(), database_name)
    filenames = list.files(input_path, pattern = "\\.csv$")
    print(filenames)
    for (fname in filenames) {
      print(fname)
      tablename = sub('.csv', '', fname)
      print(tablename)
      df = read.csv(file.path(input_path, fname), header = TRUE)
      dbWriteTable(conn, tablename, df, overwrite = TRUE, row.names = FALSE)
}
    dbDisconnect(conn)
}


load_xls = function(input_path="..//xls", database_name="summer.db"){
    # load excel spreadsheet from input_path
    conn <- dbConnect(SQLite(), database_name)
    available_sheets = c('default_constants', 'country_constants', 'default_programs', 'country_programs', 'bcg_2014', 'bcg_2015',
       'bcg_2016', 'rate_birth_2014', 'rate_birth_2015', 'life_expectancy_2014', 'life_expectancy_2015',
        'notifications_2014', 'notifications_2015', 'notifications_2016', 'outcomes_2013', 'outcomes_2015',
        'mdr_2014', 'mdr_2015', 'mdr_2016', 'laboratories_2014', 'laboratories_2015', 'laboratories_2016',
        'strategy_2014', 'strategy_2015', 'strategy_2016', 'diabetes', 'gtb_2015', 'gtb_2016', 'latent_2016',
        'tb_hiv_2016', 'spending_inputs')


     filenames = list.files(input_path, pattern = "\\.xlsx$")

     for (fname in filenames) {
       tablename = sub('.xlsx', '', fname)
       print(tablename)
       wb = loadWorkbook(file.path(input_path, fname))
       sheet_names = getSheets(wb)
       for (worksheet in sheet_names){
         if (worksheet %in% available_sheets){
           df = readWorksheet(wb, sheet = worksheet)
           dbWriteTable(conn, tablename, df, overwrite = TRUE, row.names = FALSE)
         }
       }
     }
    dbDisconnect(conn)
}

db_query = function(table_name, is_filter="", value="", column="*"){
    # method to query table_name
    conn <- dbConnect(SQLite(), "..//databases//Inputs.db")
    query = paste("Select" , column,  "from", table_name )
    if (is_filter != "" & value != ""){
        query = paste(query, "Where", is_filter,  "=", shQuote(value)  )
    }
    return(dbFetch(dbSendQuery(conn, query)))
}



db_store = function(table_name, outputs){
    # method to store outputs
    conn <- dbConnect(SQLite(), "summer.db")
    dbWriteTable(conn, table_name, outputs, overwrite = TRUE, row.names =FALSE )
    dbDisconnect(conn)
}

source('sir_runner.R')
sir_model$compartment_values
outputs  = sir_model$outputs

db_query(table_name='gtb_2015', column='year')
#conn <- dbConnect(SQLite(), "summer.db")
#bcg_2015  = dbReadTable(conn, "bcg_2015")
#res = dbSendQuery(conn, "SELECT * FROM bcg_2015 WHERE Cname = 'Bhutan'")
#bhutan_bcg = dbFetch(res)
#print(bhutan_bcg)
#dbDisconnect(conn)