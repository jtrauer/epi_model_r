library(RSQLite)
conn <- dbConnect(SQLite(), "summer.db")

filenames = list.files('xls')

for (fname in filenames) {
  tablename = str_remove(fname, '.csv')
  print(tablename)
  df = read.csv(file.path('xls', fname), header = TRUE)
  dbWriteTable(conn, tablename, df, row.names = FALSE)
}

#notification_2019 = read.csv('xls/notifications_2019.csv', header = TRUE)

#dbWriteTable(conn, "notification_2019", notification_2019, row.names = FALSE)

#head(dbReadTable(conn, "notification_2019"))


conn <- dbConnect(SQLite(), "summer.db")
source('sir_runner.R')
sir_model$compartment_values
outputs  = sir_model$outputs
dbWriteTable(conn, "model_outputs", outputs, row.names =FALSE )
dbDisconnect(conn)
