FROM node:18.16.0 as build-stage
WORKDIR /app
COPY ./ .
RUN yarn install
RUN yarn build

# FROM nginx as production-stage
# RUN mkdir /app
# COPY --from=build-stage /app/dist /app
# COPY nginx.conf /etc/nginx/nginx.conf