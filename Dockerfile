FROM alpine/git
WORKDIR /app
RUN git clone https://github.com/PankratzLab/NGS-PCA.git

FROM maven:3.6.3-jdk-8
WORKDIR /app/ngspca
COPY --from=0 /app/NGS-PCA /app
RUN mvn install

FROM openjdk:8-jre
WORKDIR /app
COPY --from=1 /app/ngspca/target/ngspca-0.01-SNAPSHOT.jar /app
ENTRYPOINT ["java", "-jar", "ngspca-0.01-SNAPSHOT.jar"]
