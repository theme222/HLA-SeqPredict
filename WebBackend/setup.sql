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
    start_time integer NOT NULL,
    duration integer NOT NULL,
    FOREIGN KEY (user_id) REFERENCES Users(id)
);

-- CODE : 0 - HAS NOT BEEN CALLED, 1 - WORKING, 2 - DONE
CREATE TABLE "User_sequences" (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    upload_time DATETIME DEFAULT CURRENT_TIMESTAMP,
    user_id INTEGER NOT NULL,
    label TEXT NOT NULL,
    igv INTEGER DEFAULT 0,
    hla_la INTEGER DEFAULT 0,
    hisat_genotype INTEGER DEFAULT 0,
    optitype INTEGER DEFAULT 0,
    FOREIGN KEY (user_id) REFERENCES Users(id)
);