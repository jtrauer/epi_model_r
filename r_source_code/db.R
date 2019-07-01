library(RSQLite)
library(XLConnect)
conn <- dbConnect(SQLite(), "summer.db")

filenames = list.files('xls', pattern = "\\.csv$")

for (fname in filenames) {
  tablename = sub('.csv', '', fname)
  print(tablename)
  df = read.csv(file.path('xls', fname), header = TRUE)
  dbWriteTable(conn, tablename, df, overwrite = TRUE, row.names = FALSE)
}

#notification_2019 = read.csv('xls/notifications_2019.csv', header = TRUE)

#dbWriteTable(conn, "notification_2019", notification_2019, row.names = FALSE)

#head(dbReadTable(conn, "notification_2019"))

available_sheets = c('default_constants', 'country_constants', 'default_programs', 'country_programs', 'bcg_2014', 'bcg_2015', 
   'bcg_2016', 'rate_birth_2014', 'rate_birth_2015', 'life_expectancy_2014', 'life_expectancy_2015',
   'notifications_2014', 'notifications_2015', 'notifications_2016', 'outcomes_2013', 'outcomes_2015',
   'mdr_2014', 'mdr_2015', 'mdr_2016', 'laboratories_2014', 'laboratories_2015', 'laboratories_2016',
   'strategy_2014', 'strategy_2015', 'strategy_2016', 'diabetes', 'gtb_2015', 'gtb_2016', 'latent_2016',
   'tb_hiv_2016', 'spending_inputs')


filenames = list.files('xls', pattern = "\\.xlsx$")

for (fname in filenames) {
  tablename = sub('.xlsx', '', fname)
  print(tablename)
  wb = loadWorkbook(file.path('xls', fname))
  sheet_names = getSheets(wb)
  for (worksheet in sheet_names){
    if (worksheet %in% available_sheets){
      df = readWorksheet(wb, sheet = worksheet)
      dbWriteTable(conn, tablename, df, overwrite = TRUE, row.names = FALSE)
    }
  }
}


conn <- dbConnect(SQLite(), "summer.db")
source('sir_runner.R')
sir_model$compartment_values
outputs  = sir_model$outputs
dbWriteTable(conn, "model_outputs", outputs, overwrite = TRUE, row.names =FALSE )
dbDisconnect(conn)

conn <- dbConnect(SQLite(), "summer.db")
bcg_2015  = dbReadTable(conn, "bcg_2015")
res = dbSendQuery(conn, "SELECT * FROM bcg_2015 WHERE Cname = 'Bhutan'")
bhutan_bcg = dbFetch(res)
print(bhutan_bcg)
dbDisconnect(conn)