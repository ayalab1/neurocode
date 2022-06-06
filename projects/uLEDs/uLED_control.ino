/* Send a TTL pulse on any of the 1-12 trigger-in channel of the OSCLite 1 for 
/ controlling the uLED, according to some specification */
// settings
int delayBetweenPulsesFrom = 20; // in ms, first value on range of for random pulses delay, for example 20
int delayBetweenPulsesTo = 40; // in ms, last value on range of for random pulses delay, for example 40
int pulseDuration = 20; // in ms, for example 20
boolean continousStimulation = true;// stimulation no dependent of the behavoiur arduino
boolean stimulateShank1 = true;
boolean stimulateShank2 = true;
boolean stimulateShank3 = true;
boolean stimulateShank4 = true; 

// Pins
const int s1l1Pin = 34;
const int s1l2Pin = 32;
const int s1l3Pin = 36;

const int s2l1Pin = 30;
const int s2l2Pin = 38;
const int s2l3Pin = 28;

const int s3l1Pin = 40;
const int s3l2Pin = 26;
const int s3l3Pin = 42;

const int s4l1Pin = 24;
const int s4l2Pin = 44;
const int s4l3Pin = 22;

const int BlockPin = 53; // From the control behaviour arduino

// state variables
int stimBlock = 0; // 0 is no stimulating
int uled_channel = 0;
int delayBetweenPulses;

void setup() {
  pinMode(LED_BUILTIN, OUTPUT);
  
  pinMode(s1l1Pin, OUTPUT);
  pinMode(s1l2Pin, OUTPUT);
  pinMode(s1l3Pin, OUTPUT);
  pinMode(s2l1Pin, OUTPUT);
  pinMode(s2l2Pin, OUTPUT);
  pinMode(s2l3Pin, OUTPUT);
  pinMode(s3l1Pin, OUTPUT);
  pinMode(s3l2Pin, OUTPUT);
  pinMode(s3l3Pin, OUTPUT);
  pinMode(s4l1Pin, OUTPUT);
  pinMode(s4l2Pin, OUTPUT);
  pinMode(s4l3Pin, OUTPUT);

  pinMode(BlockPin, INPUT);

  delayBetweenPulsesTo = delayBetweenPulsesTo + 1;
  randomSeed(analogRead(A0));   
}

void loop() {
  if (digitalRead(BlockPin) == HIGH | continousStimulation){
    digitalWrite(LED_BUILTIN, HIGH);   // turn the LED on   
    uled_channel = random(1,13); // random number from 1 to 12
    if (uled_channel == 1 & stimulateShank1 == true){
      digitalWrite(s1l1Pin, HIGH); delay(pulseDuration); digitalWrite(s1l1Pin, LOW);
      } else if (uled_channel == 2 & stimulateShank1 == true){
      digitalWrite(s1l2Pin, HIGH); delay(pulseDuration); digitalWrite(s1l2Pin, LOW);
      } else if (uled_channel == 3 & stimulateShank1 == true){
      digitalWrite(s1l3Pin, HIGH); delay(pulseDuration); digitalWrite(s1l3Pin, LOW);
      } else if (uled_channel == 4 & stimulateShank2 == true){
      digitalWrite(s2l1Pin, HIGH); delay(pulseDuration); digitalWrite(s2l1Pin, LOW);
      } else if (uled_channel == 5 & stimulateShank2 == true){
      digitalWrite(s2l2Pin, HIGH); delay(pulseDuration); digitalWrite(s2l2Pin, LOW);
      } else if (uled_channel == 6 & stimulateShank2 == true){
      digitalWrite(s2l3Pin, HIGH); delay(pulseDuration); digitalWrite(s2l3Pin, LOW);
      } else if (uled_channel == 7 & stimulateShank3 == true){
      digitalWrite(s3l1Pin, HIGH); delay(pulseDuration); digitalWrite(s3l1Pin, LOW);
      } else if (uled_channel == 8 & stimulateShank3 == true){
      digitalWrite(s3l2Pin, HIGH); delay(pulseDuration); digitalWrite(s3l2Pin, LOW);
      } else if (uled_channel == 9 & stimulateShank3 == true){
      digitalWrite(s3l3Pin, HIGH); delay(pulseDuration); digitalWrite(s3l3Pin, LOW);
      } else if (uled_channel == 10 & stimulateShank4 == true){
      digitalWrite(s4l1Pin, HIGH); delay(pulseDuration); digitalWrite(s4l1Pin, LOW);
      } else if (uled_channel == 11 & stimulateShank4 == true){
      digitalWrite(s4l2Pin, HIGH); delay(pulseDuration); digitalWrite(s4l2Pin, LOW);
      } else if (uled_channel == 12 & stimulateShank4 == true){
      digitalWrite(s4l3Pin, HIGH); delay(pulseDuration); digitalWrite(s4l3Pin, LOW);
      }
    delayBetweenPulses = random(delayBetweenPulsesFrom,delayBetweenPulsesTo); // random number from 1 to 12
    delay(delayBetweenPulses);
      
    } else {
      digitalWrite(LED_BUILTIN, LOW);   // turn the LED off   
    }
}
