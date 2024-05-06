-- FOR hla_adr.db --

CREATE TABLE "Alleles" (
    allele_id integer,
    gene varchar(255), -- A,B,C,DRB1 etc.
    allelegroup varchar(255), -- 01 02 14 etc.
    protein varchar(255),
    PRIMARY KEY (allele_id)
);

CREATE TABLE "Drugs" (
    drug_id integer,
    name varchar(255),
    PRIMARY KEY (drug_id)
);

CREATE TABLE "ADR" (
    adr_id integer,
    name varchar(255),
    PRIMARY KEY (adr_id)
);

CREATE TABLE "HLA_ADR" (
    id integer,
    allele_id integer,
    drug_id integer,
    adr_id integer,
    PRIMARY KEY (id)
    FOREIGN KEY (allele_id) REFERENCES Alleles(allele_id)
    FOREIGN KEY (drug_id) REFERENCES Drugs(drug_id)
    FOREIGN KEY (adr_id) REFERENCES ADR(adr_id)
);