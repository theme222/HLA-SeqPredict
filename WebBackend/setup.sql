-- FOR users.db --

CREATE TABLE "Users" (
    id integer PRIMARY KEY AUTOINCREMENT,
    name varchar(255),
    email varchar(255) UNIQUE,
    password_hash text NOT NULL
);

CREATE TABLE "Sessions" (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    user_id INTEGER NOT NULL,
    token TEXT NOT NULL,
    start_time INTEGER NOT NULL, -- Linux Timestamp
    duration INTEGER NOT NULL, -- Time (seconds)
    FOREIGN KEY (user_id) REFERENCES Users(id)
);

CREATE TABLE "User_sequences" (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    upload_time DATETIME DEFAULT CURRENT_TIMESTAMP,
    user_id INTEGER NOT NULL,
    label TEXT NOT NULL,
    is_paired INTEGER DEFAULT 0,  -- sqlite does not have a bool value so just using 0 and 1
    status TEXT DEFAULT "IDLE",
    igv TEXT DEFAULT "PENDING",
    hla_la TEXT DEFAULT "PENDING",
    FOREIGN KEY (user_id) REFERENCES Users(id)
);

-- Current tool status list : PENDING, RUNNING, ERROR, COMPLETED
-- Current sequence status list : IDLE, RUNNING {tool}
