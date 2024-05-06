

function test(){
    console.log('help')
    /*
    let db = new sqlite3.Database('../../resources/HLA-ADR.db', sqlite3.OPEN_READONLY, (err) => {
        if (err){
            alert('FUCK')
            console.error(err)
        }
        console.log('Connected to database')
    })

    db.all('SELECT * FROM Drugs',[],(err, rows) => {
        if (err){
            console.error(err.message)
            return
        }
        for (let row of rows){
            console.log(row)
        }
    })
    db.close()
    */
}

export {test}
