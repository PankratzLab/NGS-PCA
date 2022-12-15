FROM alpine/git
WORKDIR /app
RUN git clone https://github.com/PankratzLab/NGS-PCA.git
WORKDIR /app/NGS-PCA

FROM maven:3.6.3-openjdk-11
WORKDIR /app/ngspca
COPY --from=0 /app/NGS-PCA /app
RUN mvn install

FROM maven:3.6.3-openjdk-11
WORKDIR /app
COPY --from=1 /app/ngspca/target/ngspca-0.02-SNAPSHOT.jar /app
ENV JAVA_TOOL_OPTIONS "-XX:+UseContainerSupport -XX:MaxRAMPercentage=90.0"
ENTRYPOINT ["java", "-jar","ngspca-0.02-SNAPSHOT.jar"]
