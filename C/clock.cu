#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

struct clock {
	char* name;
	long totalTime;
	long startTime;
	struct clock* next;
};

struct clock *clocks = NULL;
long time();

void startClock(char* name) {
	struct clock *cp = clocks;

	while (cp != NULL) {
		if (strcmp(cp->name,name) == 0) {
			cp->startTime = time();
			return;
		}
		cp = cp->next;
	}
	// if you are here, no match
	cp = (struct clock*)malloc(sizeof(struct clock));
	cp->name = (char*) malloc(strlen(name)+1);
	strcpy(cp->name,name);
	cp->totalTime = 0;
	cp->startTime = time();
	cp->next = clocks;
	clocks = cp;
	return;
}

void stopClock(char* name) {
	struct clock *cp = clocks;
	while (cp && strcmp(cp->name,name)) {
		cp = cp->next;
	}
	if (cp && cp->startTime) {
		cp->totalTime += (time() - cp->startTime);
		cp->startTime = 0;
	}
}

void dump() {
	struct clock *cp = clocks;

	while (cp) {
		printf("%-20s %ld micros\n",cp->name, cp->totalTime);
		cp = cp->next;
	}
}

void printClock(char* name) {
	struct clock *cp = clocks;
	while (cp && strcmp(cp->name,name)) {
		cp = cp->next;
	}
	if (cp) {
		printf("%-20s %ld micros\n",cp->name,cp->totalTime);
	}
}

long time() {
	struct timeval tv;

	gettimeofday(&tv,NULL);
	return 1000000*(tv.tv_sec % (60*60*24*365)) + tv.tv_usec;
}

