library(RSQLite)
conn <- dbConnect(SQLite(), "summer.db")

notification_2019 = read.csv('xls/notifications_2019.csv', header = TRUE)

dbWriteTable(conn, "notification_2019", notification_2019, row.names = FALSE)

head(dbReadTable(conn, "notification_2019"))

dbDisconnect(conn)
