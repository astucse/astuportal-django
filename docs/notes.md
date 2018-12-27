# postgres


    kpostgres
    su - postgres
        create db: /usr/bin/createdb bedrock
    psql
        \c bedrock : connect to database
        CREATE USER kibru WITH PASSWORD 'supersecret';
        GRANT ALL PRIVILEGES ON DATABASE astuportal to kibru;



# pip
    pip freeze > requirements.txt
    


